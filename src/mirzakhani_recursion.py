import scipy.special
import sympy as sp
import pickle
import src.utils as utils
import scipy
import time
import logging

logger = logging.getLogger(__name__)

class WeilPetersonTable:
    def __init__(self, pickled_table):        
        self.load_table(pickled_table)

    def load_table(self, filename):
        logger.info(f"Loading table from <{filename}>.")
        with open(filename, "rb") as file:
            self.table = pickle.load(file)
        
    def save_table(self, filename):
        logger.info(f"Saving table to <{filename}>.")
        with open(filename, "wb") as file:
            pickle.dump(self.table, file)
    
    def _check_table(self, g, n):
        try:
            V = self.table[f"g={g}"][f"n={n}"]
            
            if V is None:
                found = False
            else:
                found = True
            
            logger.info(f"⋅ checking table: V_({g},{n}) - {found}")
            return V, True
        
        except KeyError:
            logger.info(f"⋅ checking table: V_({g},{n}) - {False}")
            return None, False
    
    def _add_to_table(self, g, n, V):
        key_g = f"g={g}"
        key_n = f"n={n}"
        
        if key_g not in self.table.keys():
            self.table[key_g] = {}
        self.table[key_g][key_n] = V
        
        logger.info(f"⋅ stored to table: V_({g},{n})")
        
    def initialize_table(self, filename="Weil-Peterson-base.pkl"):
        L1 = sp.Symbol("L1", positive=True)

        V_table = {
            "g=0": {"n=3": sp.Integer(1), 
                    "n=2": sp.Integer(0), 
                    "n=1": sp.Integer(0)},
            "g=1": {"n=1": L1**2 / 48 + sp.pi**2 / 12}}

        with open(filename, "wb") as f:
            pickle.dump(V_table, f)

        logger.info(f"initialized new table at <{filename}>")
        
    def display_table(self):
        logger.info("Displaying table...")
        import pandas as pd
        table = pd.DataFrame(self.table)
        table = ~table.isnull()
        table = table.reindex(sorted(table.index), axis=0)
        print()
        print(table.to_markdown(tablefmt="rounded_grid"))
        print()
        
class WeilPetersonCalculator(WeilPetersonTable):
    def __init__(self, pickled_table, exact=True):
        super().__init__(pickled_table)
        self.x = sp.Symbol("x", positive=True)      # integration variables
        self.y = sp.Symbol("y", positive=True)      # integration variables
        self.exact = exact
        
        if self.exact:
            self.factorial = sp.factorial
            self.zeta = sp.zeta
        else:
            self.factorial = scipy.special.factorial
            self.zeta = scipy.special.zeta
    
    def F(self, k, t):
        """
        Computes F_{2k-1}(t)
        """
        sum = sp.Float(0)
        k = int(k)
        for i in range(0, k+1):
            coeff = self.zeta(2*i)
            coeff *= 2**(2*i + 1) - 4
            coeff /= self.factorial(2*k-2*i)
            
            term = coeff * t**(2*k-2*i)
            sum += term.expand()
        
        sum *= self.factorial(2*k-1)
        return sum
    
    def _H(self, x, y):
        """
        returns H(x,y)
        """    
        return sp.Function("H")(x, y)
    
    def _compute_term_1(self, g, n):
        """
        requires n≧1 and g≧1, so that 2g+n-1>2 for contributions to exist
        """
        if g < 1 or n < 0:
        #if g < 1 or 1n < 0:
            return sp.Integer(0)
        
        V, found = self._check_table(g-1, n+1)

        if not found:
            V = self.calculate_V(g-1, n+1)
           
        L1 = self.L_list[0]
        x = self.x
        y = self.y
        
        logger.debug(f"({g},{n}) TERM 1")
        logger.debug(f"----------------------")
        logger.debug(f"V_({g-1},{n+1}) = {V}")
        
        V = V.subs({L1: x, self.L_list[n]: y})
        logger.debug(f"subbed: {V}")
        logger.debug(f"----------------------")
        logger.debug(f"INTEGRATING TERM 1 ({g},{n})")
        
        integrand = x * y * self._H(x + y, L1) * V
        """
        if n==0:
            V, found = self._check_table(g-1, n+1)
            integrand = x * y * self._H(x + y, sp.pi) * V.subs({L1: sp.pi})
        """
        term1  = self._double_integral(integrand, x, y)
             
        return term1
    
    def _compute_term_2(self, g, n):
        """
        Computes the third (surface-splitting) term in Mirzakhani"s recursion.
        Splits the surface into two surfaces of genus g₁ & g₂ with n₁ & n₂ boundaries, such that g₁+g₂=g and n₁-n₂= n+1.

        :param n: Number of boundaries
        :param g: Genus
        :param L_list: List of boundary length symbols
        :return: SymPy expression representing term3
        """
        if g<0: #or n<2:
            return sp.Integer(0)
        if 2*g + n  < 3:
            return sp.Integer(0)
        
        x = self.x
        y = self.y
        L1 = self.L_list[0]

        partitions = utils.all_bipartitions(self.L_list[1:n])

        integrand = sp.Integer(0)
        for g1 in range(0, g+1):
            g2 = g - g1
            for L_I, L_J in partitions:
                n1 = len(L_I)
                n2 = len(L_J)
                       
                # Sanity check the partitions: 
                if set(L_I).intersection(L_J) != set():
                    logger.warning(f"Overlapping boundaries in partition: {L_I}, {L_J}")
                if set(L_I).union(L_J) != set(self.L_list[1:n]):
                    logger.warning(f"Partition does not cover all boundaries: {L_I}, {L_J}")
                
                # Compute V₁ and V₂:
                if (2*g1 + n1 >= 2) and (2*g2 + n2 >= 2):
                    V1, found1 = self._check_table(g1, n1+1)
                    V2, found2 = self._check_table(g2, n2+1)
                    
                    if not found1:
                        V1 = self.calculate_V(g1, n1+1)
                    
                    if not found2:
                        V2 = self.calculate_V(g2, n2+1)
                    
                    V1 = V1.subs({L1: x})      
                    V1 = V1.subs({old: new for old, new in zip(self.L_list[1:n1+1], L_I)}, simultaneous=True)
                    V2 = V2.subs({L1: y})                        
                    V2 = V2.subs({old: new for old, new in zip(self.L_list[1:n2+1], L_J)}, simultaneous=True)
                        
                    integrand += x*y*self._H(x + y, L1) * V1 * V2

        term2 = self._double_integral(integrand, x, y)
        
        return term2
    
    def _compute_term_3(self, g, n):
        """
        Computes the third term in Mirzakhani's recursion.

        :param n: Number of boundaries
        :param g: Genus
        :param L_list: List of boundary length symbols
        :return: SymPy expression representing term3
        """
        if n < 2:
            return sp.Integer(0)
        
        if 2*g + n < 3:
            return sp.Integer(0)
        
        V, found = self._check_table(g, n-1)
        if not found:
            V = self.calculate_V(g, n-1)#, L_list[:n-1])
        
        x  = self.x
        L1 = self.L_list[0]
        Ln = self.L_list[n-1]
        integrand = sp.Integer(0)

        for k in range(1, n):
            Lk = self.L_list[k]
            Lk_hat = list(self.L_list[:n].copy())
            Lk_hat.pop(k)
            V_temp = V.subs({L1:x, Lk:Ln}, simultaneous=True)
            integrand += x *(self._H(x, L1 + Lk) + self._H(x, L1 - Lk)) * V_temp
        
        term3 = self._single_integral(integrand, x) 
        
        return term3
        
    def _single_integral(self, integrand, x):
        """
        Computes integral by identifying with known F_2k-1(x) functions
        """
        if integrand == 0:
            return 0
        terms = sp.Add.make_args(integrand.expand())
        result = sp.Float(0)

        k = sp.Wild("k")                   # power in x^(2k-1)
        t = sp.Wild("t")                   # argument of H
        c = sp.Wild("c", exclude=(x,))     # constant coefficient of integrand
        H = self._H(x, t)
        sub_expr = c*self.x**(2*k - 1) * H
                
        for term in terms:
            match = term.match(sub_expr)
            
            if match and match[k].is_Number:
                result += match[c]*self.F(match[k], match[t])
                
            else:
                logger.warning(f"Unmatched term in single integral: {term}")
                logger.warning(f"Matched: {match}")
                
        return result           
 
    def _double_integral(self, integrand, x, y):
        """
        Computes integral by identifying with known F_2k-1(x) functions
        """
        if integrand == 0:
            return 0
        
        terms = sp.Add.make_args(integrand.expand())
        result = sp.Float(0)
        
        a = sp.Wild("a")                        # power in x^(2a-1)
        b = sp.Wild("b")                        # power in y^(2b-1)
        t = sp.Wild("t")                        # 2nd argument of H
        c = sp.Wild("c", exclude=(x, y))        # constant coefficient integrand
        H = self._H(x + y, t)
        sub_expr = c * x**(2*a-1) * y**(2*b-1) * H
        
        for term in terms:
            match = term.match(sub_expr)

            if match and match[a].is_Number and match[b].is_Number:
                _a = int(match[a])
                _b = int(match[b])
                _t = match[t]
                _c = match[c]
                
                F_coeff  = self.factorial(2*_a - 1)
                F_coeff *= self.factorial(2*_b - 1)
                F_coeff /= self.factorial(2*_a + 2*_b - 1)            
                new = _c*F_coeff*self.F(_a + _b, _t)
                result += new
            else:
                logger.warning(f"Unmatched term in double integral: {term}")
                logger.warning(f"Matched: {match}")
            
        return result
    
    def _apply_mirzakhanis_recursion(self, g, n):#, L_list):
        """
        Computes sum A + Ad + B for Mirzakhani"s recursion
        """
        if 2 - 2*g - n > 0:
            logger.warning(f"({g},{n}) - Invalid input")
            return sp.Integer(0)
        
        else:
            term1 = self._compute_term_1(g, n)    
            term2 = self._compute_term_2(g, n)
            term3 = self._compute_term_3(g, n)
            
            logger.debug(f"({g},{n}) TERM1: {term1}")
            logger.debug(f"({g},{n}) TERM2: {term2}")
            logger.debug(f"({g},{n}) TERM3: {term3}")
            
            return term1 + term2 + term3

    def _apply_dilaton_equation(self, g, n=0):
        """
        Applies the dilaton equation to compute V_(g,n):
        (See Corollary 23 of Do's paper)
        """
        V_next, found = self._check_table(g, n+1)
        
        if not found:
            V_next = self.calculate_V(g, n+1)
        
        dilaton_lhs = V_next.diff(self.L_list[n]).subs({self.L_list[n]: 2*sp.pi*sp.I})
        
        V = dilaton_lhs / (2*sp.pi*sp.I* (2*g - 2 + n) )
        return V        
    
    def calculate_V(self, g, n):#, L_list):
        V, found = self._check_table(g, n)
        
        if not found:
            logger.info(f"⋅ computing V_({g},{n})")
            T0 = time.time()
            
            if 2*g + n <= 3:
                V = sp.Integer(0)
                
            # Apply dilaton equation if n=0
            elif n==0:
                logger.info(f". Applying Dilaton equation...")
                V = self._apply_dilaton_equation(g, n)

            # Otherwise, use Mirzakhani's recursion
            else:
                logger.info(f"⋅ Applying Mirzakhani's recursion...")
                integrand = self._apply_mirzakhanis_recursion(g, n)
                integrand = integrand.expand()
                
                T1 = time.time()
                logger.info(f"  (took  {T1-T0:.2f} s)")
                T0 = time.time()
                logger.info(f"⋅ Integrating result...")

                L1 = self.L_list[0]
                V  = sp.integrate(integrand, L1) / ( 2*L1 )

            T1 = time.time()
            logger.info(f"  (took  {T1-T0:.2f} s)")

            # Expand & add to table            
            V = V.expand()
            self._add_to_table(g, n, V)
        return V
        
    def __call__(self, g, n):
        """
        Computes Wein-Peterson volume V_(g,n) using Mirzakhani"s recursion

        :param n: no. of boundaries
        :type n: int
        :param g: genus
        :type g: int
        :returns: computed term
        :rtype: sp.Expr
        """
        T0 = time.time()
        logger.info(f"Starting recursion for V_({g},{n})")

        self.L_list = [sp.Symbol(f"L{i}", positive=True) for i in range(1, 3*g + n+1)]
        V = self.calculate_V(g, n)
        
        T1 = time.time()
        logger.info(f"Finished recursion for V_({g},{n}) - {T1-T0} s")
        
        return V
    
if __name__=="__main__":
    import time
    from datetime import datetime
    import argparse
    from pathlib import Path

    
    TODAY = datetime.now().strftime("%d-%m-%Y")
    PROJECT_ROOT = Path(__file__).parent.parent
    DATA_PATH = PROJECT_ROOT / "data"
    LOG_PATH = PROJECT_ROOT / "logs"
    print(f"{LOG_PATH / TODAY}.log")
    
    logging.info("Started.")
    logging.basicConfig(filename=f"log.log",
                        level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    
    parser = argparse.ArgumentParser(
        description="Compute Weil-Peterson volumes using Mirzakhani's recursion")
        
    tab_dir = "tables/"
    parser.add_argument("-r", "--run"       , type=bool, default=True)
    parser.add_argument("-g", "--genus"     , type=int , default=False)
    parser.add_argument("-n", "--boundaries", type=int , default=False)
    parser.add_argument("-e", "--exact"     , type=bool, default=True)
    
    parser.add_argument("-o" , "--output"      , type=str , default=False)
    parser.add_argument("-i" , "--input"       , type=str , default=f"Weil-Peterson-base.pkl")
    parser.add_argument("-io", "--input-output", type=str , default=False)
    
    parser.add_argument("-t", "--test"      , type=bool, default=False)
    parser.add_argument("-d", "--display"   , type=str , default=False)
    parser.add_argument("-v", "--verbose"   , type=bool, default=True)
    parser.add_argument("-s", "--save"      , type=bool, default=True)
    parser.add_argument("-new", "--new"     , type=str , default=False)
    
    
    args = parser.parse_args()

    if args.new:
        table = WeilPetersonTable(args.new)
        table.initialize_table(args.new)
        
    if args.display:
        table = WeilPetersonTable(args.display)
        table.display_table()

    if not args.genus and not args.boundaries:
        exit()
    
    if args.input_output:
        args.input = args.input_output
        args.output = args.input_output
        
    if args.test:
        import test_functions as test
        test.test_F()
        test.test_calculator()

    if args.run:
        calculator = WeilPetersonCalculator(f"{DATA_PATH}/{args.input}", exact=args.exact)
        
        g = args.genus
        n = args.boundaries
        V = calculator(g=g, n=n)

        if args.verbose:
            print("--------------------------------"*2)
            print(f"V_({g},{n}):")
            print("--------------------------------"*2)
            print(V.simplify())
            print("--------------------------------"*2)

        if args.save or args.output:
            if args.output:
                calculator.save_table(args.output)

            else:
                from datetime import datetime
                
                SAVE_PATH = DATA_PATH / TODAY
                
                if not SAVE_PATH.exists():
                    SAVE_PATH.mkdir(parents=True, exist_ok=True)
                
                FILENAME = datetime.now().strftime(f"table_%d-%m-%Y_%H-%M")
                calculator.save_table(f"{SAVE_PATH}/{FILENAME}.pkl")
        