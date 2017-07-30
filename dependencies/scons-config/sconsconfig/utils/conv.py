def to_iter(var):
    if hasattr(var, '__iter__') and not isinstance(var, dict):
        return var
    else:
        return [var]
