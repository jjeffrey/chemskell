import Data.List
import Data.Char
import Data.Maybe

-- ELEMENTAL DATA

elements = 
	[
		("hydrogen", "H", 1.0079),
		("helium", "He", 4.0026),
		("lithium", "Li", 6.941),
		("beryllium", "Be", 9.012),
		("boron", "B", 10.811),
		("carbon", "C", 12.011),
		("nitrogen", "N", 14.007),
		("oxygen", "O", 16.00),
		("fluorine", "F", 19.00),
		("neon", "Ne", 20.179),
		("sodium", "Na", 22.99),
		("magnesium", "Mg", 24.30),
		("aluminum", "Al", 26.98),
		("silicon", "Si", 26.98),
		("phosphorus", "P", 30.974),
		("sulfur", "S", 32.06),
		("chlorine", "Cl", 35.453),
		("argon", "Ar", 39.948),
		("potassium", "K", 39.10),
		("calcium", "Ca", 40.08),
		("scandium", "Sc", 44.96),
		("titanium", "Ti", 47.90),
		("vanadium", "V", 50.94),
		("chromium", "Cr", 52.00),
		("manganese", "Mn", 54.938),
		("iron", "Fe", 55.85),
		("cobalt", "Co", 58.93),
		("nickel", "Ni", 58.69),
		("copper", "Cu", 63.55),
		("zinc", "Zn", 65.39),
		("gallium", "Ga", 69.72),
		("germanium","Ge", 72.59),
		("arsenic", "As", 74.92),
		("selenium", "Se", 78.96),
		("bromine", "Br", 79.90),
		("krypton", "Kr", 83.80)
	]

-- ELEMENT INFO FUNCTIONS

atomicNumber :: String -> Maybe Int
atomicNumber elem = do
				elemIndex <- findIndex (\(name,symbol,massNum) -> map toLower elem == name || elem == symbol) elements
				return (elemIndex + 1)

relativeAtomicMass :: String -> Double
relativeAtomicMass elem = foldr (\(name,symbol,rAM) val-> if map toLower elem == name || elem == symbol
															then rAM
															else val) 0.0 elements

atomicWeight = relativeAtomicMass

cPeriod :: String -> Maybe Int
cPeriod elem = maybe Nothing cPeriod' $ atomicNumber elem

cPeriod' :: Int -> Maybe Int
cPeriod' aNum | aNum <= 0 = Nothing
			 | aNum <= 2 = Just 1
			 | aNum <= 10 = Just 2
			 | aNum <= 18 = Just 3
			 | aNum <= 36 = Just 4
			 | aNum <= 54 = Just 5
			 | aNum <= 86 = Just 6
			 | aNum <= 118 = Just 7
			 | aNum > 118 = Nothing


cGroup :: String -> Maybe Int
cGroup = maybe Nothing cGroup' . atomicNumber

cGroup' :: Int -> Maybe Int
cGroup' aNum | aNum <= 0 = Nothing
			 | aNum <= 1 = Just 1
			 | aNum <= 2 = Just 18
			 | aNum <= 4 = Just (aNum - 2)
			 | aNum <= 10 = Just (aNum + 8)
			 | aNum <= 12 = Just (aNum - 10)
			 | aNum <= 18 = Just aNum
			 | aNum <= 36 = Just (aNum - 18)
			 | aNum <= 54 = Just (aNum - 36)
			 | aNum <= 56 = Just (aNum - 54)
			 | aNum <= 71 = Just 3
			 | aNum <= 86 = Just (aNum - 68)
			 | aNum <= 88 = Just (aNum - 86)
			 | aNum <= 103 = Just 3
			 | aNum <= 118 = Just (aNum - 100)
			 | aNum > 118 = Nothing

cGroupCAS :: String -> Maybe String
cGroupCAS = maybe Nothing groupToCAS . maybe Nothing cGroup' . atomicNumber

groupToNames :: [String] -> Int -> Maybe String
groupToNames names num | num < 1 = Nothing
					   | num <= length names = Just $ names !! (num - 1)
					   | num > length names = Nothing

groupToCAS = groupToNames ["IA","IIA","IIIB","IVB","VB","VIB","VIB","VIIIB","VIIIB","VIIIB","IB","IIB","IIB","IIIA","IVA","VA","VIA","VIA","VIIIA"]


-- FORMULA INFO FUNCTIONS

--Gets Formula Mass From CFormula
formulaMass :: CFormula -> Double
formulaMass formula = foldr (\(symbol, num) acc -> acc + (atomicWeight symbol * fromIntegral num)) 0.0 (stitchFormula formula)

percentByMass :: CFormula -> [(CSymbol,Double)]
percentByMass formula = let totalMass = formulaMass formula in 
	map (\(symbol,num) -> (symbol, (100.0 * atomicWeight symbol * fromIntegral num) / totalMass)) (stitchFormula formula)

-- Example Formulas:

-- Multiple Level Example
-- CComp [(CUnit [("Ba",1)],1),( CComp[(CUnit [("O",1),("H",1)], 1)],2)]

-- Single Level Example
-- CComp [(CUnit [("Ba",1)],1),(CUnit [("O",1),("H",1)],2)]


-- RAW FORMULA FUNCTIONS (currently cannot handle charges)

type CSymbol = String
data CFormula = CUnit [(CSymbol, Int)] | CComp [(CFormula, Int)] deriving Eq

instance Show CFormula where
	show (CUnit unit) = foldr (\(symb, num) str -> symb ++ (if num > 1 then show num else "") ++ str) "" unit
	show (CComp comp) = foldr (\(formula, num) str -> (if num > 1 then "(" ++ show formula ++ ")" ++ show num else show formula) ++ str) "" comp

data CEquation = CEquation ([(Int, CFormula)], [(Int, CFormula)]) deriving (Eq, Show)

--Converts CComp into CUnit (if CUnit, returns CUnit)
simplifyFormula :: CFormula -> CFormula
simplifyFormula = CUnit . stitchFormula

--Properly extracts data from CFormula
stitchFormula :: CFormula -> [(CSymbol, Int)]
stitchFormula (CComp comp) = foldr (\(formula, num) stitched -> stitchFormula formula ++ stitched) [] ((\(CComp comp) -> comp) $ flattenFormula $ CComp comp)
stitchFormula (CUnit unit) = unit

--Makes all CComp coefficients 1 and scales CUnits appropriately
flattenFormula :: CFormula -> CFormula
flattenFormula (CComp comp) = CComp $ map (\(formula, num) -> (flattenFormula (scaleFormula formula num),1)) comp
flattenFormula (CUnit unit) = CUnit unit

-- Shallow scaling (does not go to deepest level in CComp)
scaleFormula :: CFormula -> Int -> CFormula
scaleFormula (CComp comp) scale = CComp $ map (\(formula,num) -> (formula, num * scale)) comp
scaleFormula (CUnit unit) scale = CUnit $ map (\(symb,num) -> (symb, num * scale)) unit


-- VALIDATION FUNCTIONS

--Pattern-matching version
validateCSymbol :: CSymbol -> Bool
validateCSymbol symbol = (isUpper . head) symbol && length symbol < 4 && all isLower (tail symbol)

--uses whitelisting of symbols based on elements
isElementSymbol :: CSymbol -> Bool
isElementSymbol symbol = any (\(name, symb, massNum) -> symbol == symb) elements

--uses whitelisting based on elements
isElement :: String -> Bool
isElement str = any (\(name, symbol, massNum) -> name == map toLower str || symbol == str) elements


-- PARSING FUNCTIONS

formulaFromStructure :: [String] -> CFormula
formulaFromStructure structure = undefined


-- finds distinct symbols, coefficients, and parentheses in a string (Assumes string is written properly)
getStructure :: String -> [String]
getStructure = 
	foldr (\a acc-> 
		if isUpper a && acc /= [] && all isLower (head acc)
			then if drop 2 (head acc) /= ""
				then (a : take 2 (head acc)) : drop 2 (head acc) : tail acc
				else (a : take 2 (head acc)) : tail acc
			else if (isDigit a && acc/= [] && all isDigit (head acc)) || (isLower a && acc/= [] && all isLower (head acc))
				then (a : head acc) : tail acc
				else [a] : acc
	) []

--Function Needed to Traverse from structure to CFormula
matchingOpenParenIndex :: [String] -> Int
matchingOpenParenIndex structure = 
	1 + snd (foldr (\a (count, indx) -> 
			if count < 1 
				then (if a == ")" 
					then count - 1 
					else if a == "(" 
						then count + 1 
						else count, indx - 1)
				else (2, indx))
		(0, length structure - 1) structure) :: Int