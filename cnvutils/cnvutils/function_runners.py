import inspect
import itertools
import multiprocessing

from .load_data import load_event_metadata

def _runner(
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

    print(f"Running {func.__name__}...\n")# with {kwargs}...\n")

    return func(**kwargs)

def multi_runner(
    func,
    sources,
    levels,
    chromosomes_events,
    more_dicts=[],
    threads=None,
):

    args = _get_multi_runner_args(
        func=func,
        sources=sources,
        levels=levels,
        chromosomes_events=chromosomes_events,
        more_dicts=more_dicts,
    )

    if threads is None:
        with multiprocessing.Pool() as pool:
            results = pool.starmap(_runner, args)
    else:
        with multiprocessing.Pool(threads) as pool:
            results = pool.starmap(_runner, args)

    return results

def _get_multi_runner_args(
    func,
    sources,
    levels,
    chromosomes_events,
    more_dicts,
):

    # Save each combination of args we want to run
    args = []
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

                        # For any additional parameters, just do all combinations of them
                        if more_dicts:
                            names = [more_dict["name"] for more_dict in more_dicts]
                            vals = [more_dict["vals"] for more_dict in more_dicts]
                            for combo in itertools.product(*vals):
                                more_args = dict(zip(names, combo))
                                args.append((func, meta, more_args))

                        else:
                            more_args = {}
                            args.append((func, meta, more_args))
    return args
