class Borg:
    _shared_state = {}

    def __new__(cls):
        instance = super().__new__(cls)
        instance.__dict__ = cls._shared_state
        return instance