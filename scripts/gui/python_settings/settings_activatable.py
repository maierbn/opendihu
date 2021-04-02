class Activatable:
    activated = True
    parent = None
    def __init__(self):
        pass

    def activate_recursive(self):
        self.activated = True
        if self.parent:
            self.parent.activate_recursive()
