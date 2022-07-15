// $ANTLR 3.4 /ews/QSPR-Core/src/qspr/oscript/OScript.g 2012-02-21 18:31:35
/* Copyright (C) 2022 BIGCHEM GmbH <info@bigchem.de>
 *
 * Contact: info@bigchem.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License (AGPL)
 * as published by the Free Software Foundation; either version 3.0
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the Affero GNU General Public License for more details.
 *
 * You should have received a copy of the Affero GNU Lesser General Public License
 * along with this program; If not, see <https://www.gnu.org/licenses/>. 
 */

package qspr.oscript; 

import org.antlr.runtime.BitSet;
import org.antlr.runtime.EarlyExitException;
import org.antlr.runtime.MismatchedSetException;
import org.antlr.runtime.NoViableAltException;
import org.antlr.runtime.Parser;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.RecognizerSharedState;
import org.antlr.runtime.Token;
import org.antlr.runtime.TokenStream;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class OScriptParser extends Parser {
	public static final String[] tokenNames = new String[] {
			"<invalid>", "<EOR>", "<DOWN>", "<UP>", "ARROW", "COMMENT", "IDENTIFIER", "IDENTIFIER_CHAR", "NEWLINE", "NUMBER", "PREDICATE", "WHITESPACE", "'*'", "'['", "']'"
	};

	public static final int EOF=-1;
	public static final int T__12=12;
	public static final int T__13=13;
	public static final int T__14=14;
	public static final int ARROW=4;
	public static final int COMMENT=5;
	public static final int IDENTIFIER=6;
	public static final int IDENTIFIER_CHAR=7;
	public static final int NEWLINE=8;
	public static final int NUMBER=9;
	public static final int PREDICATE=10;
	public static final int WHITESPACE=11;

	// delegates
	public Parser[] getDelegates() {
		return new Parser[] {};
	}

	// delegators


	public OScriptParser(TokenStream input) {
		this(input, new RecognizerSharedState());
	}
	public OScriptParser(TokenStream input, RecognizerSharedState state) {
		super(input, state);
	}

	public String[] getTokenNames() { return OScriptParser.tokenNames; }
	public String getGrammarFileName() { return "/ews/QSPR-Core/src/qspr/oscript/OScript.g"; }



	// $ANTLR start "start"
	// /ews/QSPR-Core/src/qspr/oscript/OScript.g:6:1: start returns [TransformationRule transformationRule] : ( ( rule | COMMENT ) ( WHITESPACE )* ( NEWLINE )* ( WHITESPACE )* )+ EOF ;
	public final TransformationRule start() throws RecognitionException {
		TransformationRule transformationRule = null;


		AtomicRule rule1 =null;


		try {
			// /ews/QSPR-Core/src/qspr/oscript/OScript.g:6:54: ( ( ( rule | COMMENT ) ( WHITESPACE )* ( NEWLINE )* ( WHITESPACE )* )+ EOF )
			// /ews/QSPR-Core/src/qspr/oscript/OScript.g:7:2: ( ( rule | COMMENT ) ( WHITESPACE )* ( NEWLINE )* ( WHITESPACE )* )+ EOF
			{
				transformationRule = new TransformationRule();

				// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:2: ( ( rule | COMMENT ) ( WHITESPACE )* ( NEWLINE )* ( WHITESPACE )* )+
				int cnt5=0;
				loop5:
					do {
						int alt5=2;
						int LA5_0 = input.LA(1);

						if ( ((LA5_0 >= COMMENT && LA5_0 <= IDENTIFIER)||LA5_0==12) ) {
							alt5=1;
						}


						switch (alt5) {
						case 1 :
							// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:3: ( rule | COMMENT ) ( WHITESPACE )* ( NEWLINE )* ( WHITESPACE )*
						{
							// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:3: ( rule | COMMENT )
							int alt1=2;
							int LA1_0 = input.LA(1);

							if ( (LA1_0==IDENTIFIER||LA1_0==12) ) {
								alt1=1;
							}
							else if ( (LA1_0==COMMENT) ) {
								alt1=2;
							}
							else {
								NoViableAltException nvae =
										new NoViableAltException("", 1, 0, input);

								throw nvae;

							}
							switch (alt1) {
							case 1 :
								// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:4: rule
							{
								pushFollow(FOLLOW_rule_in_start34);
								rule1=rule();

								state._fsp--;


								transformationRule.addRule(rule1);

							}
							break;
							case 2 :
								// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:58: COMMENT
							{
								match(input,COMMENT,FOLLOW_COMMENT_in_start38); 

							}
							break;

							}


							// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:67: ( WHITESPACE )*
							loop2:
								do {
									int alt2=2;
									int LA2_0 = input.LA(1);

									if ( (LA2_0==WHITESPACE) ) {
										alt2=1;
									}


									switch (alt2) {
									case 1 :
										// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:67: WHITESPACE
									{
										match(input,WHITESPACE,FOLLOW_WHITESPACE_in_start41); 

									}
									break;

									default :
										break loop2;
									}
								} while (true);


							// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:79: ( NEWLINE )*
							loop3:
								do {
									int alt3=2;
									int LA3_0 = input.LA(1);

									if ( (LA3_0==NEWLINE) ) {
										alt3=1;
									}


									switch (alt3) {
									case 1 :
										// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:79: NEWLINE
									{
										match(input,NEWLINE,FOLLOW_NEWLINE_in_start44); 

									}
									break;

									default :
										break loop3;
									}
								} while (true);


							// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:88: ( WHITESPACE )*
							loop4:
								do {
									int alt4=2;
									int LA4_0 = input.LA(1);

									if ( (LA4_0==WHITESPACE) ) {
										alt4=1;
									}


									switch (alt4) {
									case 1 :
										// /ews/QSPR-Core/src/qspr/oscript/OScript.g:8:88: WHITESPACE
									{
										match(input,WHITESPACE,FOLLOW_WHITESPACE_in_start47); 

									}
									break;

									default :
										break loop4;
									}
								} while (true);


						}
						break;

						default :
							if ( cnt5 >= 1 ) break loop5;
							EarlyExitException eee =
									new EarlyExitException(5, input);
							throw eee;
						}
						cnt5++;
					} while (true);


				match(input,EOF,FOLLOW_EOF_in_start54); 

			}

		}
		catch (RecognitionException re) {
			reportError(re);
			recover(input,re);
		}

		finally {
			// do for sure before leaving
		}
		return transformationRule;
	}
	// $ANTLR end "start"



	// $ANTLR start "rule"
	// /ews/QSPR-Core/src/qspr/oscript/OScript.g:11:1: rule returns [AtomicRule atomicRule] :ruleText= ( (l= specifier |i= IDENTIFIER ) ARROW (r= specifier |i= IDENTIFIER ) ) ;
	public final AtomicRule rule() throws RecognitionException {
		AtomicRule atomicRule = null;


		Token i=null;
		Token ruleText=null;
		Specifier l =null;

		Specifier r =null;


		try {
			// /ews/QSPR-Core/src/qspr/oscript/OScript.g:11:38: (ruleText= ( (l= specifier |i= IDENTIFIER ) ARROW (r= specifier |i= IDENTIFIER ) ) )
			// /ews/QSPR-Core/src/qspr/oscript/OScript.g:12:2: ruleText= ( (l= specifier |i= IDENTIFIER ) ARROW (r= specifier |i= IDENTIFIER ) )
			{
				atomicRule = new AtomicRule();

				// /ews/QSPR-Core/src/qspr/oscript/OScript.g:13:11: ( (l= specifier |i= IDENTIFIER ) ARROW (r= specifier |i= IDENTIFIER ) )
				// /ews/QSPR-Core/src/qspr/oscript/OScript.g:14:3: (l= specifier |i= IDENTIFIER ) ARROW (r= specifier |i= IDENTIFIER )
				{
					// /ews/QSPR-Core/src/qspr/oscript/OScript.g:14:3: (l= specifier |i= IDENTIFIER )
					int alt6=2;
					int LA6_0 = input.LA(1);

					if ( (LA6_0==IDENTIFIER) ) {
						int LA6_1 = input.LA(2);

						if ( (LA6_1==13) ) {
							alt6=1;
						}
						else if ( (LA6_1==ARROW) ) {
							alt6=2;
						}
						else {
							NoViableAltException nvae =
									new NoViableAltException("", 6, 1, input);

							throw nvae;

						}
					}
					else if ( (LA6_0==12) ) {
						alt6=1;
					}
					else {
						NoViableAltException nvae =
								new NoViableAltException("", 6, 0, input);

						throw nvae;

					}
					switch (alt6) {
					case 1 :
						// /ews/QSPR-Core/src/qspr/oscript/OScript.g:14:4: l= specifier
					{
						pushFollow(FOLLOW_specifier_in_rule81);
						l=specifier();

						state._fsp--;


						atomicRule.left = l;

					}
					break;
					case 2 :
						// /ews/QSPR-Core/src/qspr/oscript/OScript.g:14:53: i= IDENTIFIER
					{
						i=(Token)match(input,IDENTIFIER,FOLLOW_IDENTIFIER_in_rule89); 

						atomicRule.left = new Specifier((i!=null?i.getText():null));

					}
					break;

					}


					match(input,ARROW,FOLLOW_ARROW_in_rule97); 

					// /ews/QSPR-Core/src/qspr/oscript/OScript.g:16:3: (r= specifier |i= IDENTIFIER )
					int alt7=2;
					int LA7_0 = input.LA(1);

					if ( (LA7_0==IDENTIFIER) ) {
						int LA7_1 = input.LA(2);

						if ( (LA7_1==13) ) {
							alt7=1;
						}
						else if ( (LA7_1==EOF||(LA7_1 >= COMMENT && LA7_1 <= IDENTIFIER)||LA7_1==NEWLINE||(LA7_1 >= WHITESPACE && LA7_1 <= 12)) ) {
							alt7=2;
						}
						else {
							NoViableAltException nvae =
									new NoViableAltException("", 7, 1, input);

							throw nvae;

						}
					}
					else if ( (LA7_0==12) ) {
						alt7=1;
					}
					else {
						NoViableAltException nvae =
								new NoViableAltException("", 7, 0, input);

						throw nvae;

					}
					switch (alt7) {
					case 1 :
						// /ews/QSPR-Core/src/qspr/oscript/OScript.g:16:4: r= specifier
					{
						pushFollow(FOLLOW_specifier_in_rule104);
						r=specifier();

						state._fsp--;


						atomicRule.right = r;

					}
					break;
					case 2 :
						// /ews/QSPR-Core/src/qspr/oscript/OScript.g:16:54: i= IDENTIFIER
					{
						i=(Token)match(input,IDENTIFIER,FOLLOW_IDENTIFIER_in_rule112); 

						atomicRule.right = new Specifier((i!=null?i.getText():null));

					}
					break;

					}


				}


				atomicRule.rawRule = (ruleText!=null?ruleText.getText():null);

			}

		}
		catch (RecognitionException re) {
			reportError(re);
			recover(input,re);
		}

		finally {
			// do for sure before leaving
		}
		return atomicRule;
	}
	// $ANTLR end "rule"



	// $ANTLR start "specifier"
	// /ews/QSPR-Core/src/qspr/oscript/OScript.g:21:1: specifier returns [Specifier specifier] : (property= ( IDENTIFIER | '*' ) '[' condition= IDENTIFIER PREDICATE value= ( IDENTIFIER | NUMBER ) ( WHITESPACE unit= IDENTIFIER )? ']' ) ;
	public final Specifier specifier() throws RecognitionException {
		Specifier specifier = null;


		Token property=null;
		Token condition=null;
		Token value=null;
		Token unit=null;
		Token PREDICATE2=null;

		try {
			// /ews/QSPR-Core/src/qspr/oscript/OScript.g:21:42: ( (property= ( IDENTIFIER | '*' ) '[' condition= IDENTIFIER PREDICATE value= ( IDENTIFIER | NUMBER ) ( WHITESPACE unit= IDENTIFIER )? ']' ) )
			// /ews/QSPR-Core/src/qspr/oscript/OScript.g:22:2: (property= ( IDENTIFIER | '*' ) '[' condition= IDENTIFIER PREDICATE value= ( IDENTIFIER | NUMBER ) ( WHITESPACE unit= IDENTIFIER )? ']' )
			{
				// /ews/QSPR-Core/src/qspr/oscript/OScript.g:22:2: (property= ( IDENTIFIER | '*' ) '[' condition= IDENTIFIER PREDICATE value= ( IDENTIFIER | NUMBER ) ( WHITESPACE unit= IDENTIFIER )? ']' )
				// /ews/QSPR-Core/src/qspr/oscript/OScript.g:22:3: property= ( IDENTIFIER | '*' ) '[' condition= IDENTIFIER PREDICATE value= ( IDENTIFIER | NUMBER ) ( WHITESPACE unit= IDENTIFIER )? ']'
				{
					property=(Token)input.LT(1);

					if ( input.LA(1)==IDENTIFIER||input.LA(1)==12 ) {
						input.consume();
						state.errorRecovery=false;
					}
					else {
						MismatchedSetException mse = new MismatchedSetException(null,input);
						throw mse;
					}


					match(input,13,FOLLOW_13_in_specifier150); 

					condition=(Token)match(input,IDENTIFIER,FOLLOW_IDENTIFIER_in_specifier154); 

					PREDICATE2=(Token)match(input,PREDICATE,FOLLOW_PREDICATE_in_specifier156); 

					value=(Token)input.LT(1);

					if ( input.LA(1)==IDENTIFIER||input.LA(1)==NUMBER ) {
						input.consume();
						state.errorRecovery=false;
					}
					else {
						MismatchedSetException mse = new MismatchedSetException(null,input);
						throw mse;
					}


					// /ews/QSPR-Core/src/qspr/oscript/OScript.g:22:91: ( WHITESPACE unit= IDENTIFIER )?
					int alt8=2;
					int LA8_0 = input.LA(1);

					if ( (LA8_0==WHITESPACE) ) {
						alt8=1;
					}
					switch (alt8) {
					case 1 :
						// /ews/QSPR-Core/src/qspr/oscript/OScript.g:22:92: WHITESPACE unit= IDENTIFIER
					{
						match(input,WHITESPACE,FOLLOW_WHITESPACE_in_specifier167); 

						unit=(Token)match(input,IDENTIFIER,FOLLOW_IDENTIFIER_in_specifier171); 

					}
					break;

					}


					match(input,14,FOLLOW_14_in_specifier175); 

				}


				specifier = new Specifier((property!=null?property.getText():null), (condition!=null?condition.getText():null), (PREDICATE2!=null?PREDICATE2.getText():null), (value!=null?value.getText():null), (unit!=null?unit.getText():null));

			}

		}
		catch (RecognitionException re) {
			reportError(re);
			recover(input,re);
		}

		finally {
			// do for sure before leaving
		}
		return specifier;
	}
	// $ANTLR end "specifier"

	// Delegated rules




	public static final BitSet FOLLOW_rule_in_start34 = new BitSet(new long[]{0x0000000000001960L});
	public static final BitSet FOLLOW_COMMENT_in_start38 = new BitSet(new long[]{0x0000000000001960L});
	public static final BitSet FOLLOW_WHITESPACE_in_start41 = new BitSet(new long[]{0x0000000000001960L});
	public static final BitSet FOLLOW_NEWLINE_in_start44 = new BitSet(new long[]{0x0000000000001960L});
	public static final BitSet FOLLOW_WHITESPACE_in_start47 = new BitSet(new long[]{0x0000000000001860L});
	public static final BitSet FOLLOW_EOF_in_start54 = new BitSet(new long[]{0x0000000000000002L});
	public static final BitSet FOLLOW_specifier_in_rule81 = new BitSet(new long[]{0x0000000000000010L});
	public static final BitSet FOLLOW_IDENTIFIER_in_rule89 = new BitSet(new long[]{0x0000000000000010L});
	public static final BitSet FOLLOW_ARROW_in_rule97 = new BitSet(new long[]{0x0000000000001040L});
	public static final BitSet FOLLOW_specifier_in_rule104 = new BitSet(new long[]{0x0000000000000002L});
	public static final BitSet FOLLOW_IDENTIFIER_in_rule112 = new BitSet(new long[]{0x0000000000000002L});
	public static final BitSet FOLLOW_set_in_specifier143 = new BitSet(new long[]{0x0000000000002000L});
	public static final BitSet FOLLOW_13_in_specifier150 = new BitSet(new long[]{0x0000000000000040L});
	public static final BitSet FOLLOW_IDENTIFIER_in_specifier154 = new BitSet(new long[]{0x0000000000000400L});
	public static final BitSet FOLLOW_PREDICATE_in_specifier156 = new BitSet(new long[]{0x0000000000000240L});
	public static final BitSet FOLLOW_set_in_specifier160 = new BitSet(new long[]{0x0000000000004800L});
	public static final BitSet FOLLOW_WHITESPACE_in_specifier167 = new BitSet(new long[]{0x0000000000000040L});
	public static final BitSet FOLLOW_IDENTIFIER_in_specifier171 = new BitSet(new long[]{0x0000000000004000L});
	public static final BitSet FOLLOW_14_in_specifier175 = new BitSet(new long[]{0x0000000000000002L});

}