import sympy as sp
from src import utils 
from src.mirzakhani_recursion import WeilPetersonCalculator
from pathlib import Path
import logging

def test_F():
    TEST_PATH = Path(__file__).parent
    tester = WeilPetersonCalculator(pickled_table = TEST_PATH / "test_table.pkl")
    t = sp.Symbol("t")

    # Expected results for k=1, k=2 and k=3:
    F_1 = t**2 / 2 + 2*sp.pi**2/3
    F_3 = t**4/4 + 2*sp.pi**2*t**2 + 28*sp.pi**4/15
    F_5 = t**6/6 + 10*sp.pi**2*t**4/3 + 56*sp.pi**4*t**2/3 + 992*sp.pi**6/63

    expected = [F_1, F_3, F_5]    
    computed = [tester.F(k, t) for k in range(1, 4)]
    
    for i in range(3):
        msg = f"Computed F_{2*i+1}(t) did not match expected:"
        msg += f"\nExpected: {expected[i]}"
        msg += f"\nComputed: {computed[i]}"
           
        assert abs(expected[i] - computed[i])<1e-14, msg

def test_V06(instance):
    L    = [sp.Symbol(f"L{i}", positive=True) for i in range(1,7)]
    m3   = utils.m([3], L)
    m21  = utils.m([2,1], L)
    m111 = utils.m([1,1,1], L)   
    m2   = utils.m([2], L)
    m11  = utils.m([1,1], L)
    m1   = utils.m([1], L)

    expected = (m3/48 + m21*3/16 + m111*3/4 + m2*3*sp.pi**2/2 + m11*6*sp.pi**2 +
                m1*26*sp.pi**4 + 244*sp.pi**6/3).expand()
    computed = instance(g=0,n=6).as_expr()
    
    assert computed.equals(expected), "Test failed for (g,n) = (0,6)"
    
def test_V05(instance):
    L   = [sp.Symbol(f"L{i}", positive=True) for i in range(1,6)]
    m2  = utils.m([2], L)
    m11 = utils.m([1,1], L)
    m1  = utils.m([1], L)
    
    expected = (m2/8 + m11/2 + m1*3*sp.pi**2 + 10*sp.pi**4).expand()
    computed = instance(g=0,n=5).as_expr()
    
    assert computed.equals(expected), "Test failed for (g,n) = (0,5)"

def test_V04(instance):
    L  = [sp.Symbol(f"L{i}", positive=True) for i in range(1,5)]
    m1 = utils.m([1], L)
    
    expected = (m1/2 + 2*sp.pi**2).expand()
    computed = instance(g=0, n=4).as_expr()
    
    assert computed.equals(expected), "Test failed for (g,n) = (0,4)"

def test_V13(instance):
    L    = [sp.Symbol(f"L{i}", positive=True) for i in range(1,4)]
    m3   = utils.m([3], L)
    m21  = utils.m([2,1], L)
    m111 = utils.m([1,1,1], L)
    m2   = utils.m([2], L)
    m11  = utils.m([1,1], L)
    m1   = utils.m([1], L)
    
    expected = (m3/1152 + m21/192 + m111/96 + m2*sp.pi**2/24 + 
               m11*sp.pi**2/8 + m1*13*sp.pi**4/24 + 14*sp.pi**6/9).expand()
    computed = instance(g=1,n=3).as_expr()
    
    assert computed.equals(expected), "Test failed for (g,n) = (1,3)"

def test_V12(instance):
    L   = [sp.Symbol(f"L{i}", positive=True) for i in range(1,3)]
    m2  = utils.m([2], L)
    m11 = utils.m([1,1], L)
    m1  = utils.m([1], L)
    
    expected = (m2/192 + m11/96 + m1*sp.pi**2/12 +sp.pi**4 / 4).expand()
    computed = instance(g=1, n=2).as_expr()
    
    #print(expected)
    #print(computed)
    assert computed.equals(expected), "Test failed for (g,n) = (1,2)"

def test_V14(instance):
    L = [sp.Symbol(f"L{i}", positive=True) for i in range(1,4+1)]
    m1    = utils.m([1], L)
    m11   = utils.m([1,1], L)
    m111  = utils.m([1,1,1], L)
    m1111 = utils.m([1,1,1,1], L)
    m2    = utils.m([2], L)
    m21   = utils.m([2,1], L)
    m211  = utils.m([2,1,1], L)
    m22   = utils.m([2,2], L)
    m3    = utils.m([3], L)
    m31   = utils.m([3,1], L)
    m4    = utils.m([4], L)
    
    expected = (m4/9216 + m31/768 + m22/384 + m211/128 + m1111/64 + m3*7*sp.pi**2/576 
                + m21*sp.pi**2/12 + m111*sp.pi**2/4 + m2*41*sp.pi**4/96 + m11*17*sp.pi**4/12 + 
                m1*187*sp.pi**6/36 + 529*sp.pi**8/36).expand()
    computed = instance(g=1, n=4).as_expr()
    
    assert computed.equals(expected), "Test failed for (g,n) = (1,4)"

def test_calculator():
    TEST_PATH = Path(__file__).parent
    tester = WeilPetersonCalculator(pickled_table = TEST_PATH / "test_table.pkl")
    
    test_V06(tester)
    test_V05(tester)
    test_V04(tester)
    test_V12(tester)
    test_V13(tester)
    test_V14(tester)
    
    logging.info("All tests passed!")
    
if __name__ == "__main__":
    test_F()
    test_calculator()
