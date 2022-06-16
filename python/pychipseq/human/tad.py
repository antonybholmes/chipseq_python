from .. import tad
from . import human


class GencodeTADAnnotation(tad.TADAnnotation):
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
        
    def __init__(self, bin_size=100):
        super().__init__(human.TAD_FILE, bin_size=bin_size)