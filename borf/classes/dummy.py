class Dummy:
    def __init__(self):
        self.text = "text"

    def __str__(self):
        return self.text
    __repr__ = __str__