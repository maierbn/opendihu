import copy
from root_node import RootNode
from helpers import Error, Info

# holds a list of different RootNodes and can switch between them to undo and redo actions


class UndoStack:
    def __init__(self):
        self.stack = []
        self.current_index = -1

    def get_current_root(self):
        return self.stack[self.current_index]

    def duplicate_current_state(self):
        # deepcopy current root
        self.remove_future()
        self.stack.append(copy.deepcopy(self.get_current_root()))
        # swap the new copy with the current root (so current root is at the end of the stack)
        self.stack[self.current_index], self.stack[self.current_index +
                                                   1] = self.stack[self.current_index + 1], self.stack[self.current_index]
        self.current_index = self.current_index + 1

    def undo(self):
        if self.current_index > 0:
            self.current_index = self.current_index - 1
            return Info('undo successful')
        else:
            return Error('cannot undo')

    def redo(self):
        if len(self.stack) - 1 > self.current_index:
            self.current_index = self.current_index + 1
            return Info('redo successful')
        else:
            return Error('cannot redo')

    def add(self, node):
        self.remove_future()
        self.stack.append(node)
        self.current_index = self.current_index + 1

    def remove_future(self):
        # pop everything newer than the current root
        self.stack = self.stack[:self.current_index + 1]
