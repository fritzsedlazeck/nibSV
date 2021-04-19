# vim: sw=4 ts=4 sts=4 tw=0 et:
from nibpkg/welcome import nil
import unittest

suite "welcome":
    test "home":
        assert 1 == 1
        assert welcome.getWelcomeMessage() == "Hello, World!"
