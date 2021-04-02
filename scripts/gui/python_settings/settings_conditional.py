# if-else block inside a SettingsDictEntry.value or a SettingsListEntry.value
class SettingsConditional():
    def __init__(self):
        self.condition = None
        self.if_block = None
        self.else_block = None

    def repr(self, depth, hide_placeholders=False):
        depth = depth - 1
        if isinstance(self.if_block, str):
            value1 = self.if_block
        else:
            value1 = self.if_block.repr(depth + 1, hide_placeholders=hide_placeholders)
        if isinstance(self.else_block, str):
            value2 = self.else_block
        else:
            value2 = self.else_block.repr(depth + 1, hide_placeholders=hide_placeholders)
        return value1 + ' if ' + self.condition + ' else ' + value2
