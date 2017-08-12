
class MyClass:
    """A simple example class"""
    i = 12345
    
    def f(self):
        return 'hello world'

# import mytest
# x=mytest.MyClass.i
# x=mytest.MyClass.f()

class Dog:

    kind = 'canine'         # class variable shared by all instances

    def __init__(self, name):
        self.name = name    # instance variable unique to each instance

# x=mytest.Dog('Fido')
# x.name
