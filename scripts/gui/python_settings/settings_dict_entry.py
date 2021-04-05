from python_settings.settings_activatable import Activatable

# normal entry in a SettingsDict
class SettingsDictEntry(Activatable):
    def __init__(self, key=None, value=None, comment=None, doc_link=None):
        if isinstance(key, str) and not key[0] == '"':
            self.key = '"' + key + '"'
        else:
            self.key = key
        self.value = value
        self.comments = []
        if comment:
            self.comments.append('# ' + comment)
        self.doc_link = None
        self.is_unknown = False
        self.parent = None
        if doc_link:
            self.doc_link = doc_link
