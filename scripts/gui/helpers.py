import sys

# use printe() instead of print() to print errors to stderr instead of stdout


def printe(message):
    print('Error: ' + message, file=sys.stderr)


class Message:
    def __init__(self, message):
        self.message = message
        if not hasattr(self, 'prefix'):
            self.prefix = 'Message'
        if not hasattr(self, 'color'):
            self.color = None
        # self.print()

    def __repr__(self):
        return str(self.prefix + ': ' + self.message)

    def print(self):
        print(self.prefix + ': ' + self.message, file=sys.stderr)


class Error(Message):
    def __init__(self, message):
        self.prefix = 'Error'
        self.color = 'red'
        super().__init__(message)


class Info(Message):
    def __init__(self, message):
        self.prefix = 'Info'
        self.color = 'green'
        super().__init__(message)


class Warning(Message):
    def __init__(self, message):
        self.prefix = 'Warning'
        self.color = 'orange'
        super().__init__(message)

# helper function to indent a multiline-string by a given indentation


def indent(lines, indentation):
    return indentation + lines.replace('\n', '\n' + indentation)
