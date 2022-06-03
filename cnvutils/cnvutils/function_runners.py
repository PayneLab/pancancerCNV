import inspect
import multiprocessing

from .load_data import load_event_metadata

def _runner(
    lock,
    func,
    meta,
    more={},
):
    kwargs = {}
    for param in inspect.signature(func).parameters.values():
        
        name = param.name
        
        if name in meta.keys():
            val = meta[name]
        elif name in more.keys():
            val = more[name]
        elif param.default is not param.empty:
            val = param.default
        else:
            raise ValueError(f"Must provide value for param '{name}'")
            
        kwargs[name] = val

    lock.acquire()
    try:
        print(f"Running {func.__name__} with {kwargs}...\n\n")
    finally:
        lock.release()

    return #func(**kwargs)

def multi_runner(
    func,
    sources,
    levels,
    chromosomes_events,
    more_dicts=[],
    _more_args={},
):

    with multiprocessing.Pool() as pool:

        lock = multiprocessing.Lock()
        args = _get_multi_runner_args(
            lock=lock,
            func=func,
            sources=sources,
            levels=levels,
            chromosomes_events=chromosomes_events,
            more_dicts=more_dicts,
            more_args=_more_args,
        )

        results = pool.starmap(_runner, args)

    return results

def _get_multi_runner_args(
    lock,
    func,
    sources,
    levels,
    chromosomes_events,
    more_dicts,
    more_args,
):
    args = []

    if more_dicts: # We need to handle additional parameter permutations via recursion
        param = more_dicts[0]
        for val in param["vals"]:
            more_args[param["name"]] = val

            args.extend(_get_multi_runner_args(
                lock=lock,
                func=func,
                sources=sources,
                levels=levels,
                chromosomes_events=chromosomes_events,
                more_dicts=more_dicts[1:],
                more_args=more_args,
            ))

    else: # We have all the parameters we need, time to run the function

        # Save each combination of args we want to run
        for source in sources:
            
            # Handle different level options for CPTAC/GISTIC
            if source == "cptac":
                source_levels = [None]
            elif source == "gistic":
                source_levels = levels
                
            for level in source_levels:
                for chromosome, arms in chromosomes_events.items():
                    for arm, types in arms.items():
                        for event_type in types:
                            
                            meta = load_event_metadata(
                                source=source,
                                chromosome=chromosome,
                                arm=arm,
                                gain_or_loss=event_type,
                                level=level,
                            )

                            args.append((lock, func, meta, more_args))
    return args
