def frozen(set):
    """
    Raise an error when trying to set an undeclared name, or when calling
    from a method other than Frozen.__init__ or the __init__ method of
    a class derived from Frozen
    """
    def set_attr(self, name, value):
        import sys
        if hasattr(self, name):
            # If attribute already exists, simply set it
            set(self, name, value)
            return
        elif sys._getframe(1).f_code.co_name is '__init__':
            # Allow __setattr__ calls in __init__ calls of proper object types
            for k, v in sys._getframe(1).f_locals.items():
                if k == "self" and isinstance(v, self.__class__):
                    set(self, name, value)
                    return
        raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr


class FrozenClass(object):
    """
    Subclasses of Frozen are frozen, i.e. it is impossible to add
    new attributes to them and their instances.
    """
    __setattr__ = frozen(object.__setattr__)

    class __metaclass__(type):
        __setattr__ = frozen(type.__setattr__)

# An alternative to the solution above. If this option is chosen,
# call 'self._freeze()' at the end of __init__() for that class.
# class FrozenClass(object):
#     __isfrozen = False

#     def __setattr__(self, key, value):
#         if self.__isfrozen and not hasattr(self, key):
#             raise TypeError("%r is a frozen class" % (self,))
#         object.__setattr__(self, key, value)

#     def _freeze(self):
#         self.__isfrozen = True
