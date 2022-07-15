grammar UnitConversion;

@header { package qspr.util.unitconversion; }
@lexer::header { package qspr.util.unitconversion; }

@members {
    public double varvalue;  // Value to-be-converted.
    public double varmw;     // Molecular weight of compound. 
    
    public static double log(double base, double arg) { return Math.log(arg) / Math.log(base); }  // Logarithm with arbitrary base.
}

/*  ************
 *  *  Parser  *
 *  ************/

// Start symbol. Recognizes a single formula expression.
start returns [double value]
    : sumexpr { $value = $sumexpr.value; } EOF ;

// Recognizes a sum of one or more summands. +,- have lowest precedence.
sumexpr returns [double value]
    : v=mulexpr { $value = $v.value; } 
          (
            ( PLUS  w=mulexpr { $value += $w.value; } ) 
          | ( MINUS w=mulexpr { $value -= $w.value; } ) 
          )* ;

// Recognizes a product of one or more factors. *,/ have second-lowest precedence.
mulexpr returns [double value]
    : v=unaryexpr { $value = $v.value; } 
          (
            ( MUL  w=unaryexpr { $value *= $w.value; } ) 
          | ( DIV  w=unaryexpr { $value /= $w.value; } ) 
          )* ;  

// Recognizes sign of unary expression. Unary signs have third-lowest precedence.
unaryexpr	returns [double value]
    : 
          ( MINUS primexpr { $value = -$primexpr.value; } )
        | ( primexpr { $value = $primexpr.value; } )
        ;

// Recognizes primary expressions (number literals, named constants, function evaluations, bracketed expressions).
primexpr returns [double value]
    :     ( number   { $value = $number.value;  } )
        | ( var      { $value = $var.value;     } )
        | ( constant { $value = $constant.value;} )
        | ( funexpr  { $value = $funexpr.value; } )
        | ( parexpr  { $value = $parexpr.value; } ) 
        ;  

// Named numerical constants.
constant	returns [double value]
    :	    ( E  { $value = 2.7182818284590452; } )
	      | ( PI { $value = 3.1415926535897932; } )
	      | ( NA { $value = 6.02214179e23;      } )
	      ;

funexpr returns [double value]
    :     ( COS  parexpr { $value = Math.cos  ($parexpr.value); } )
	      | ( SIN  parexpr { $value = Math.sin  ($parexpr.value); } )
	      | ( SQRT parexpr { $value = Math.sqrt ($parexpr.value); } )
	      | ( POW  LPAR v=sumexpr COMMA w=sumexpr RPAR { $value = Math.pow($v.value,$w.value); } )
	      | ( EXP  parexpr { $value = Math.exp  ($parexpr.value); } )
	      | ( LOG  LPAR v=sumexpr COMMA w=sumexpr RPAR { $value = log($v.value,$w.value); } )
	      | ( LOGD parexpr { $value = Math.log10($parexpr.value); } )
	      | ( LOGN parexpr { $value = Math.log  ($parexpr.value); } )
	      | ( LOGB parexpr { $value = log(2.,$parexpr.value);     } )
	      ;

parexpr	 returns [double value]
    :	LPAR sumexpr { $value = $sumexpr.value; } RPAR ;

// Recognizes the variable placeholder and returns its value.
var  returns [double value]
    :    ( VARVALUE { $value = this.varvalue; } ) 
       | ( VARMW    { $value = this.varmw;    } )
       ;

// Recognizes a number and returns its value.
number  returns [double value]
    :	NUMBER { $value = Double.parseDouble($NUMBER.text); } ;

/*  ***********
 *  *  Lexer  *
 *  ***********/

// Whitespace.

WHITESPACE : ( '\t' | ' ' | '\r' | '\n'| '\u000C' )+ { $channel = HIDDEN; } ;

// Numerical constants.

NUMBER 	: INT | FLOAT ;

fragment INT : DIGIT+ ;  // Opening zeros, +0 and -0 are allowed.
fragment FLOAT : INT '.' DIGIT* ;

fragment DIGIT	: '0'..'9' ;

// Basic arithmetic.

PLUS    : '+'     ;  // Addition.
MINUS   : '-'     ;  // Subtraction.
MUL     : '*'     ;  // Multiplication.
DIV     : '/'     ;  // Division.

// Higher mathematical functions.

COMMA   : ','     ;  // Comma separates arguments of functions.

COS     : 'cos'   ;  // Cosine.
SIN     : 'sin'   ;  // Sine.
SQRT    : 'sqrt'  ;  // Square root.    
POW     : 'pow'   ;  // Power.
EXP     : 'exp'   ;  // Exponentiation to base E.
LOG     : 'log'   ;         // Logarithm.
LOGD    : 'lg' | 'log10' ;  // Decadic logarithm.
LOGN    : 'ln' | 'logn'  ;  // Natural logarithm.
LOGB    : 'log2'  ;         // Binary logarithm.

// Named constants.

E       : 'E'     ;  // Eulers constant.
PI      : 'PI'    ;  // Pi.
NA      : 'NA'    ;  // Avogadro constant.

// Grouping of expressions.

LPAR : '(' ;  // Left parenthesis.
RPAR : ')' ;  // Right parenthesis.

// Special.
    
VARVALUE : '#'      ;  // The converted value.
VARMW    : '$MW'    ;  // The molecular weight of the compound.
