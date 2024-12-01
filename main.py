from datetime import datetime
import argparse
import logging
from pathlib import Path
from src.mirzakhani_recursion import WeilPetersonCalculator, WeilPetersonTable

TODAY = datetime.now().strftime("%d-%m-%Y")
PROJECT_ROOT = Path(__file__).parent
DATA_PATH = PROJECT_ROOT / "data"
LOG_PATH = PROJECT_ROOT / "logs"

logging.basicConfig(filename=f"{LOG_PATH / TODAY}.log",
                    level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logging.info("Session started. ")


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

if args.test:
    import test.test_functions as test
    logging.info("Running unit tests.")
    test.test_F()
    test.test_calculator()

if args.new:
    table = WeilPetersonTable(args.new)
    table.initialize_table(args.new)
    
if args.display:
    table = WeilPetersonTable(args.display)
    table.display_table()

if not args.genus and not args.boundaries:
    logging.warning("No (g,n) specified. Exiting.")
    logging.info("Session finished.\n    ___________\n")
    exit()

if args.input_output:
    args.input = args.input_output
    args.output = args.input_output
    
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
    
logging.info("Session finished.\n    ___________\n")