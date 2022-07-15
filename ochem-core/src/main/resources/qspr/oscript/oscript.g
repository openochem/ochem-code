grammar OScript;

@header { package qspr.oscript; }
@lexer::header { package qspr.oscript; }

start returns [TransformationRule transformationRule]: 
	{$transformationRule = new TransformationRule();}
	((rule {$transformationRule.addRule($rule.atomicRule);}|COMMENT) WHITESPACE* NEWLINE* WHITESPACE*)+ 
	EOF;
	
rule returns [AtomicRule atomicRule] :
	{$atomicRule = new AtomicRule();} 
	ruleText=(
		(l=specifier {$atomicRule.left = $l.specifier;} | i=IDENTIFIER {$atomicRule.left = new Specifier($i.text);}) 
		ARROW
		(r=specifier {$atomicRule.right = $r.specifier;} | i=IDENTIFIER {$atomicRule.right = new Specifier($i.text);})
	)
	{$atomicRule.rawRule = $ruleText.text;} 
	; 

specifier returns [Specifier specifier] 	: 
	(property=(IDENTIFIER|'*')  '[' condition=IDENTIFIER PREDICATE value=(IDENTIFIER|NUMBER) (WHITESPACE unit=IDENTIFIER)? ']') {$specifier = new Specifier($property.text, $condition.text, $PREDICATE.text, $value.text, $unit.text);};
	

IDENTIFIER	: (('a'..'z'|'A'..'Z'|'_') IDENTIFIER_CHAR*) | ('"' ('a'..'z'|'A'..'Z'|'_') (IDENTIFIER_CHAR|' '|')'|'(')* '"');
NEWLINE		: '\n' | '\r' ;
ARROW		:	WHITESPACE* '->' WHITESPACE*;
PREDICATE	:	WHITESPACE* ('='|'<'|'>')  WHITESPACE*;
NUMBER	:	'0'..'9'+ | '0'..'9'+ '.' '0'..'9'+;
COMMENT	:	'/''/' ~NEWLINE* NEWLINE* ;
WHITESPACE 	:'\t'|' ';
IDENTIFIER_CHAR 
	:	'a'..'z'|'A'..'Z'|'_'|'0'..'9'|'%'|'-';
