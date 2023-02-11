class RingError(ValueError):
    """ raise this when there is a probelm with rings """
    def __init__(self, msg: object) -> object:
        self.msg = msg
        super(RingError, self).__init__(msg)
