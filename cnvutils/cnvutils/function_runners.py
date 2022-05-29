import inspect

from .load_data import load_event_metadata

def runner(
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

    print(f"Running {func.__name__} with {kwargs}...")
    return func(**kwargs)

def multi_runner(
    func,
    sources,
    levels,
    chromosomes_events,
    more_dicts=[],
    _more_args={},
):
    results = []

    # Handle additional parameters
    if more_dicts:
        param = more_dicts[0]
        for val in param["vals"]:
            _more_args[param["name"]] = val

            results.extend(multi_runner(
                func=func,
                sources=sources,
                levels=levels,
                chromosomes_events=chromosomes_events,
                more_dicts=more_dicts[1:],
                _more_args=_more_args,
            ))

    # Run with standard parameters
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

                        results.append(runner(
                            func=func,
                            meta=meta,
                            more=_more_args,
                        ))

                        print()
    return results
