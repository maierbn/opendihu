import os

def get_prefix():
    import sconsconfig as cfg
    return os.path.dirname(cfg.__file__)

def get_data_prefix():
    return os.path.join(get_prefix(), 'data')
