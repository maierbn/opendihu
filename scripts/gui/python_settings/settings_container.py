import re
import typing
import sys

from tokenize import tokenize, untokenize, NUMBER, STRING, NAME, OP
from helpers import Error
import token
from io import BytesIO

from python_settings.settings_child_placeholder import *
from python_settings.settings_list_entry import *
from python_settings.settings_dict_entry import *
from python_settings.settings_activatable import *
from python_settings.settings_empty_line import *
from python_settings.settings_conditional import *
from python_settings.settings_choice import *

# this class is the parent of SettingsDict and SettingsList
class SettingsContainer(list):
    # replaces the first found SettingsChildPlaceholder with the entries of child_dict
    def replaceChildPlaceholder(self, child_dict):
        for i in range(len(self)):
            entry = self[i]
            if isinstance(entry, SettingsChildPlaceholder):
                self.pop(i)
                while len(child_dict) > 0:
                    self.insert(i, child_dict.pop())
                return True
            elif (isinstance(entry, SettingsListEntry) or isinstance(entry, SettingsDictEntry)) and isinstance(entry.value, SettingsContainer):
                if entry.value.replaceChildPlaceholder(child_dict):
                    return True

    # counts the SettingsChildPlaceholder that are direct childs of this SettingsContainer
    def count_child_placeholders(self):
        count = 0
        for entry in self:
            if isinstance(entry, SettingsChildPlaceholder):
                count = count + 1
        return count


# represents a python-settings-dict
class SettingsDict(SettingsContainer, Activatable):
    # init an empty SettingsDict or parse a settings-string to a SettingsDict
    # you can also give this a list with entries
    def __init__(self, settings=None):
        if settings == None:
            return
        elif isinstance(settings, typing.List):
            for entry in settings:
                self.append(entry)
            return

        # remove outer braces
        settings = settings[1:][:-1]

        # split settings into tokens using tokenize from the stdlib
        tokens = tokenize(BytesIO(settings.encode('utf-8')).readline)

        stack = []
        stack.append(self)
        mode_stack = []
        mode_stack.append("dict_key")
        nested_counter = 0

        token_buffer = []
        append_comment = False
        token_type_last = None
        for t in tokens:
            token_value = t.string
            token_type = token.tok_name[t.exact_type]
            #print()
            #try:
            #   print(stack[0])
            #except: pass
            #print(type(stack[-1]))
            #print(mode_stack)
            #print(token_value)
            #print(token_type + token_value)
            if token_type == 'NEWLINE' or token_type == 'NL':
                # in the edgecase where there is no comma after a value -> store the value
                if nested_counter == 0 and len(token_buffer) > 0:
                    append_comment = True
                    if isinstance(stack[-1], SettingsDict):
                        if len(token_buffer) > 0:
                            stack[-1][-1].value = tokens_to_string(
                                token_buffer)
                            token_buffer = []
                        mode_stack.pop()
                    else:
                        if len(token_buffer) > 0:
                            list_entry = SettingsListEntry()
                            list_entry.value = tokens_to_string(token_buffer)
                            stack[-1].append(list_entry)
                            token_buffer = []
                # don't append comments to SettingsDictEntry or SettingsListEntry after newline
                append_comment = False
                # handle empty lines
                if token_type_last == 'NEWLINE' or token_type_last == 'NL':
                    stack[-1].append(SettingsEmptyLine())
            # handle comments
            elif token_type == 'COMMENT':
                if append_comment:
                    # append the comment to the last list entry
                    stack[-1][-1].comments.append(token_value)
                else:
                    reg_child_placeholder = re.compile('### CHILD [0-9] ###')
                    match_child_placeholder = reg_child_placeholder.match(
                        token_value)
                    # reg_link_placeholder = re.compile('### [A-Z]+ ###')
                    #match_link_placeholder = reg_link_placeholder.match(token_value)
                    if match_child_placeholder:
                        c = SettingsChildPlaceholder(token_value[10:-4])
                    # elif match_link_placeholder:
                    #    c = SettingsLinkPlaceholder(token_value[4:-4])
                    else:
                        # add floating comment to the current SettingsDict or SettingsList
                        c = SettingsComment()
                        c.comment = token_value
                    stack[-1].append(c)
            # handle curly braces '{}'
            elif token_type == 'LBRACE':
                if mode_stack[-1] == "list_comprehension" or mode_stack[-1] == 'conditional':
                    nested_counter = nested_counter + 1
                elif nested_counter == 0:
                    mode_stack.append("dict_key")
                    dict = SettingsDict()
                    if isinstance(stack[-1], SettingsList):
                        stack[-1].append(SettingsListEntry())
                    if isinstance(mode_stack[-1], SettingsConditional):
                        stack[-1].else_block = list
                    else:
                        stack[-1][-1].value = dict
                    stack.append(dict)
                    append_comment = False
            elif token_type == 'RBRACE':
                if mode_stack[-1] == "list_comprehension" or mode_stack[-1] == 'conditional':
                    nested_counter = nested_counter - 1
                elif nested_counter == 0:
                    if len(token_buffer) > 0:
                        stack[-1][-1].value = tokens_to_string(token_buffer)
                        token_buffer = []
                    # pop 2 times because of key+value on mode_stack
                    mode_stack.pop()
                    # but if there was a trailing COMMA before, value is already popped, so we only have to pop once
                    if mode_stack[-1] == 'dict_key':
                        mode_stack.pop()
                    stack.pop()
                    append_comment = True
            # handle square brackets '[]'
            elif token_type == 'LSQB':
                if mode_stack[-1] == "list_comprehension" or mode_stack[-1] == 'conditional':
                    nested_counter = nested_counter + 1
                elif nested_counter == 0:
                    mode_stack.append("list")
                    # detected child list
                    list = SettingsList()
                    if isinstance(stack[-1], SettingsList):
                        stack[-1].append(SettingsListEntry())
                    if isinstance(stack[-1], SettingsConditional):
                        stack[-1].else_block = list
                    else:
                        stack[-1][-1].value = list
                    stack.append(list)
                    append_comment = False
            elif token_type == 'RSQB':
                if mode_stack[-1] == "list_comprehension":
                    if nested_counter == 1:
                        #if token_buffer:
                        stack[-1].list_comprehension = tokens_to_string(
                        token_buffer)
                        token_buffer = []
                        nested_counter = 0
                        stack.pop()
                        mode_stack.pop()
                        mode_stack.pop()
                    else:
                        nested_counter = nested_counter - 1
                elif mode_stack[-1] == 'conditional':
                    nested_counter = nested_counter - 1
                elif nested_counter == 0:
                    if len(token_buffer) > 0:
                        if isinstance(stack[-1], SettingsList):
                            list_entry = SettingsListEntry()
                            list_entry.value = tokens_to_string(token_buffer)
                            stack[-1].append(list_entry)
                        token_buffer = []
                    mode_stack.pop()
                    stack.pop()
                    append_comment = True
                    if isinstance(stack[-1], SettingsConditional):
                        # also pop the SettingsConditional, we are done with that
                        stack.pop()

            # handle dictionary keys
            elif mode_stack[-1] == "dict_key":
                if token_type == 'STRING' or token_type == 'NUMBER':
                    # we got a new key (keys are always STRING or NUMBER)
                    stack[-1].append(SettingsDictEntry())
                    append_comment = True
                    stack[-1][-1].key = token_value
                elif token_type == 'COLON':
                    mode_stack.append("dict_value")

            # handle dictionary values and list entries
            elif mode_stack[-1] == "dict_value" or mode_stack[-1] == "list" or mode_stack[-1] == "list_comprehension" or mode_stack[-1] == "conditional":
                # handle comma ',' (only if we are not in nested braces)
                if nested_counter == 0 and token_type == 'COMMA':
                    append_comment = True
                    if isinstance(stack[-1], SettingsDict):
                        if len(token_buffer) > 0:
                            stack[-1][-1].value = tokens_to_string(
                                token_buffer)
                            token_buffer = []
                        mode_stack.pop()
                    else:
                        if len(token_buffer) > 0:
                            list_entry = SettingsListEntry()
                            list_entry.value = tokens_to_string(token_buffer)
                            stack[-1].append(list_entry)
                            token_buffer = []
                elif nested_counter == 0 and token_type == 'NAME' and token_value == 'if':
                    nested_counter = nested_counter + 1
                    mode_stack.append('conditional')
                    # get the value of the last entry we added and replace it with a SettingsConditional containing the old value
                    first_condition_value = stack[-1][-1].value
                    stack[-1][-1].value = SettingsConditional()
                    stack[-1][-1].value.if_block = first_condition_value

                    stack.append(stack[-1][-1].value)
                elif mode_stack[-1] == 'conditional' and nested_counter == 1 and token_type == 'NAME' and token_value == 'else':
                    nested_counter = 0
                    mode_stack.pop()
                    stack[-1].condition = tokens_to_string(token_buffer)
                    token_buffer = []
                else:
                    # handle list-comprehensions like [ ... for i in range 10]
                    if nested_counter == 0 and token_type == 'NAME' and token_value == 'for':
                        nested_counter = nested_counter + 1
                        mode_stack.append('list_comprehension')

                    # handle other tokens
                    # if not already continued, token_value must be part of the value
                    token_buffer.append(t)

                    # handle parentheses '()'
                    if token_type == 'LPAR':
                        nested_counter = nested_counter + 1
                    if token_type == 'RPAR':
                        nested_counter = nested_counter - 1

            token_type_last = token_type

    def __repr__(self):
        return self.repr(0)

    def repr(self, depth, hide_placeholders=False):
        if len(self) == 0:
            return '{}'
        indentation = '  '
        r = ''
        for i in range(len(self)):
            entrie = self[i]
            # don't print inserted defaults entries, that are not activated
            try:
                if entrie.activated == False:
                    continue
            except: pass
            if isinstance(entrie, SettingsDictEntry):
                comments = ''
                for comment in entrie.comments:
                    comments = comments + ' ' + comment
                if isinstance(entrie.value, str):
                    value = entrie.value
                else:
                    #print(entrie.key)
                    #print(entrie.comments)
                    #print(entrie.value)
                    value = entrie.value.repr(depth + 1, hide_placeholders=hide_placeholders)
                optional_comma = ','
                if i == len(self) - 1:
                    optional_comma = ''
                entrie_r = indentation * \
                    (depth + 1) + entrie.key + ' : ' + \
                    value + optional_comma + comments
            elif hide_placeholders and isinstance(entrie, SettingsChildPlaceholder):
                continue
            elif isinstance(entrie, SettingsChoice) or isinstance(entrie, SettingsMesh) or isinstance(entrie, SettingsSolver):
                # repr everything in SettingsChoice SettingsMesh and SettingsSolver (only for testing, this should not happen normally)
                entrie_r = ''
                if isinstance(entrie, SettingsChoice):
                    es = entrie.defaults + entrie.alternatives
                else:
                    es = entrie
                    # resolve choices in meshes and solvers:
                    for e in es:
                        if isinstance(e, SettingsChoice):
                            es.extend(e.defaults + e.alternatives)
                            es.remove(e)
                for e in es:
                    if isinstance(e.value, str):
                        value = e.value
                    else:
                        value = e.value.repr(depth + 1, hide_placeholders=hide_placeholders)
                    entrie_r = entrie_r + indentation * (depth + 1) + e.key + ' : ' + value + ',\n'
                entrie_r = entrie_r[:-1]
            elif isinstance(entrie, SettingsComment):
                # SettingsChildPlaceholder gets handled here if hide_placeholders==False
                entrie_r = indentation * (depth + 1) + entrie.comment
            elif isinstance(entrie, SettingsEmptyLine):
                entrie_r = ''
            r = r + '\n' + entrie_r
        return '{' + r + '\n' + indentation * depth + '}'

    def has_key(self, key):
        conditionals_resolved = self.__get_resolved_Conditionals()
        return any(isinstance(entry, SettingsDictEntry) and key == entry.key for entry in self + conditionals_resolved)

    def get_value(self, key):
        conditionals_resolved = self.__get_resolved_Conditionals()
        for entry in self + conditionals_resolved:
            if isinstance(entry, SettingsDictEntry) and entry.key == key:
                return entry.value
        return

    def get_entry(self, key):
        conditionals_resolved = self.__get_resolved_Conditionals()
        for entry in self + conditionals_resolved:
            if isinstance(entry, SettingsDictEntry) and entry.key == key:
                return entry
        return

    def __get_resolved_Conditionals(self):
        conditionals_resolved = []
        for entry in self:
            if isinstance(entry, SettingsConditional):
                # this only resolves stuff in if_block not in else_block
                for e in entry.if_block:
                    conditionals_resolved.append(e)
        return conditionals_resolved

# represents a list stored in a SettingsDictEntry.value or a SettingsListEntry.value
class SettingsList(SettingsContainer, Activatable):
    def __init__(self, entries=None):
        self.list_comprehension = None
        if entries:
            if isinstance(entries, str):
                # wrap list-str in dict and parse it
                l = SettingsDict('{ "a" : ' + entries + ' }')[0].value
                self.list_comprehension = l.list_comprehension
                for e in l:
                    self.append(e)
            elif isinstance(entries, list):
                for entry in entries:
                    # if we only pass strings, we don't have to encapsualte them manually
                    if isinstance(entry, str):
                        entry = SettingsListEntry(entry)
                    self.append(entry)

    def __repr__(self):
        return self.repr(0)

    def repr(self, depth, hide_placeholders=False):
        if len(self) == 0:
            return '[]'
        indentation = '  '
        r = ''
        for i in range(len(self)):
            entrie = self[i]
            if isinstance(entrie, SettingsListEntry):
                comments = ''
                for comment in entrie.comments:
                    comments = comments + ' ' + comment
                if isinstance(entrie.value, str):
                    value = entrie.value
                else:
                    value = entrie.value.repr(depth + 1, hide_placeholders=hide_placeholders)
                optional_comma = ','
                comprehension = ''
                if i == len(self) - 1:
                    optional_comma = ''
                    if self.list_comprehension:
                        comprehension = ' ' + self.list_comprehension
                entrie_r = indentation * \
                    (depth + 1) + value + comprehension + \
                    optional_comma + comments
            elif isinstance(entrie, SettingsComment):
                entrie_r = indentation * (depth + 1) + entrie.comment
            elif isinstance(entrie, SettingsEmptyLine):
                entrie_r = ''
            r = r + '\n' + entrie_r
        return '[' + r + '\n' + indentation * depth + ']'

    def get_first_SettingsListEntry(self):
        for entry in self:
            if isinstance(entry, SettingsListEntry):
                return entry
        return

    # get the i'th SettingsListEntry
    def get_settings_list_entry(self, i):
        j = 0
        for entry in self:
            if isinstance(entry, SettingsListEntry):
                if i == j:
                    return entry
                j = j + 1
        return None


class SettingsMesh(SettingsDict):
    def __init__(self, options):
        for entry in options:
            self.append(entry)
        self.name_key = '"meshName"'
        self.name_prefix = 'mesh'
        self.global_key = '"Meshes"'

class SettingsSolver(SettingsDict):
    def __init__(self, options):
        for entry in options:
            self.append(entry)
        self.name_key = '"solverName"'
        self.name_prefix = 'solver'
        self.global_key = '"Solvers"'

# helper function wrapping pythons untokenize-function to improve readability of the returned string
def tokens_to_string(tokens):
    #print(tokens)
    return untokenize(tokens).splitlines()[-1].strip()
