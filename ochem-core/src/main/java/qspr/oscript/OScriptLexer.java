// $ANTLR 3.4 /ews/QSPR-Core/src/qspr/oscript/OScript.g 2012-02-21 18:31:36
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

import org.antlr.runtime.BaseRecognizer;
import org.antlr.runtime.CharStream;
import org.antlr.runtime.DFA;
import org.antlr.runtime.EarlyExitException;
import org.antlr.runtime.Lexer;
import org.antlr.runtime.MismatchedSetException;
import org.antlr.runtime.NoViableAltException;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.RecognizerSharedState;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class OScriptLexer extends Lexer {
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
    // delegators
    public Lexer[] getDelegates() {
        return new Lexer[] {};
    }

    public OScriptLexer() {} 
    public OScriptLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public OScriptLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);
    }
    public String getGrammarFileName() { return "/ews/QSPR-Core/src/qspr/oscript/OScript.g"; }

    // $ANTLR start "T__12"
    public final void mT__12() throws RecognitionException {
        try {
            int _type = T__12;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:4:7: ( '*' )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:4:9: '*'
            {
            match('*'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__12"

    // $ANTLR start "T__13"
    public final void mT__13() throws RecognitionException {
        try {
            int _type = T__13;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:5:7: ( '[' )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:5:9: '['
            {
            match('['); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__13"

    // $ANTLR start "T__14"
    public final void mT__14() throws RecognitionException {
        try {
            int _type = T__14;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:6:7: ( ']' )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:6:9: ']'
            {
            match(']'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__14"

    // $ANTLR start "IDENTIFIER"
    public final void mIDENTIFIER() throws RecognitionException {
        try {
            int _type = IDENTIFIER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:12: ( ( ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR )* ) | ( '\"' ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR | ' ' | ')' | '(' )* '\"' ) )
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( ((LA3_0 >= 'A' && LA3_0 <= 'Z')||LA3_0=='_'||(LA3_0 >= 'a' && LA3_0 <= 'z')) ) {
                alt3=1;
            }
            else if ( (LA3_0=='\"') ) {
                alt3=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 3, 0, input);

                throw nvae;

            }
            switch (alt3) {
                case 1 :
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:14: ( ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR )* )
                    {
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:14: ( ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR )* )
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:15: ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR )*
                    {
                    if ( (input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
                        input.consume();
                    }
                    else {
                        MismatchedSetException mse = new MismatchedSetException(null,input);
                        recover(mse);
                        throw mse;
                    }


                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:39: ( IDENTIFIER_CHAR )*
                    loop1:
                    do {
                        int alt1=2;
                        int LA1_0 = input.LA(1);

                        if ( (LA1_0=='%'||LA1_0=='-'||(LA1_0 >= '0' && LA1_0 <= '9')||(LA1_0 >= 'A' && LA1_0 <= 'Z')||LA1_0=='_'||(LA1_0 >= 'a' && LA1_0 <= 'z')) ) {
                            alt1=1;
                        }


                        switch (alt1) {
                    	case 1 :
                    	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
                    	    {
                    	    if ( input.LA(1)=='%'||input.LA(1)=='-'||(input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
                    	        input.consume();
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;
                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop1;
                        }
                    } while (true);


                    }


                    }
                    break;
                case 2 :
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:59: ( '\"' ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR | ' ' | ')' | '(' )* '\"' )
                    {
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:59: ( '\"' ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR | ' ' | ')' | '(' )* '\"' )
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:60: '\"' ( 'a' .. 'z' | 'A' .. 'Z' | '_' ) ( IDENTIFIER_CHAR | ' ' | ')' | '(' )* '\"'
                    {
                    match('\"'); 

                    if ( (input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
                        input.consume();
                    }
                    else {
                        MismatchedSetException mse = new MismatchedSetException(null,input);
                        recover(mse);
                        throw mse;
                    }


                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:25:88: ( IDENTIFIER_CHAR | ' ' | ')' | '(' )*
                    loop2:
                    do {
                        int alt2=2;
                        int LA2_0 = input.LA(1);

                        if ( (LA2_0==' '||LA2_0=='%'||(LA2_0 >= '(' && LA2_0 <= ')')||LA2_0=='-'||(LA2_0 >= '0' && LA2_0 <= '9')||(LA2_0 >= 'A' && LA2_0 <= 'Z')||LA2_0=='_'||(LA2_0 >= 'a' && LA2_0 <= 'z')) ) {
                            alt2=1;
                        }


                        switch (alt2) {
                    	case 1 :
                    	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
                    	    {
                    	    if ( input.LA(1)==' '||input.LA(1)=='%'||(input.LA(1) >= '(' && input.LA(1) <= ')')||input.LA(1)=='-'||(input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
                    	        input.consume();
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;
                    	    }


                    	    }
                    	    break;

                    	default :
                    	    break loop2;
                        }
                    } while (true);


                    match('\"'); 

                    }


                    }
                    break;

            }
            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "IDENTIFIER"

    // $ANTLR start "NEWLINE"
    public final void mNEWLINE() throws RecognitionException {
        try {
            int _type = NEWLINE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:26:10: ( '\\n' | '\\r' )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            {
            if ( input.LA(1)=='\n'||input.LA(1)=='\r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NEWLINE"

    // $ANTLR start "ARROW"
    public final void mARROW() throws RecognitionException {
        try {
            int _type = ARROW;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:27:8: ( ( WHITESPACE )* '->' ( WHITESPACE )* )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:27:10: ( WHITESPACE )* '->' ( WHITESPACE )*
            {
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:27:10: ( WHITESPACE )*
            loop4:
            do {
                int alt4=2;
                int LA4_0 = input.LA(1);

                if ( (LA4_0=='\t'||LA4_0==' ') ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            	    {
            	    if ( input.LA(1)=='\t'||input.LA(1)==' ' ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop4;
                }
            } while (true);


            match("->"); 



            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:27:27: ( WHITESPACE )*
            loop5:
            do {
                int alt5=2;
                int LA5_0 = input.LA(1);

                if ( (LA5_0=='\t'||LA5_0==' ') ) {
                    alt5=1;
                }


                switch (alt5) {
            	case 1 :
            	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            	    {
            	    if ( input.LA(1)=='\t'||input.LA(1)==' ' ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop5;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ARROW"

    // $ANTLR start "PREDICATE"
    public final void mPREDICATE() throws RecognitionException {
        try {
            int _type = PREDICATE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:28:11: ( ( WHITESPACE )* ( '=' | '<' | '>' ) ( WHITESPACE )* )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:28:13: ( WHITESPACE )* ( '=' | '<' | '>' ) ( WHITESPACE )*
            {
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:28:13: ( WHITESPACE )*
            loop6:
            do {
                int alt6=2;
                int LA6_0 = input.LA(1);

                if ( (LA6_0=='\t'||LA6_0==' ') ) {
                    alt6=1;
                }


                switch (alt6) {
            	case 1 :
            	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            	    {
            	    if ( input.LA(1)=='\t'||input.LA(1)==' ' ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop6;
                }
            } while (true);


            if ( (input.LA(1) >= '<' && input.LA(1) <= '>') ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:28:40: ( WHITESPACE )*
            loop7:
            do {
                int alt7=2;
                int LA7_0 = input.LA(1);

                if ( (LA7_0=='\t'||LA7_0==' ') ) {
                    alt7=1;
                }


                switch (alt7) {
            	case 1 :
            	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            	    {
            	    if ( input.LA(1)=='\t'||input.LA(1)==' ' ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop7;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "PREDICATE"

    // $ANTLR start "NUMBER"
    public final void mNUMBER() throws RecognitionException {
        try {
            int _type = NUMBER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:29:8: ( ( '0' .. '9' )+ | ( '0' .. '9' )+ '.' ( '0' .. '9' )+ )
            int alt11=2;
            alt11 = dfa11.predict(input);
            switch (alt11) {
                case 1 :
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:29:10: ( '0' .. '9' )+
                    {
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:29:10: ( '0' .. '9' )+
                    int cnt8=0;
                    loop8:
                    do {
                        int alt8=2;
                        int LA8_0 = input.LA(1);

                        if ( ((LA8_0 >= '0' && LA8_0 <= '9')) ) {
                            alt8=1;
                        }


                        switch (alt8) {
                    	case 1 :
                    	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
                    	        input.consume();
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;
                    	    }


                    	    }
                    	    break;

                    	default :
                    	    if ( cnt8 >= 1 ) break loop8;
                                EarlyExitException eee =
                                    new EarlyExitException(8, input);
                                throw eee;
                        }
                        cnt8++;
                    } while (true);


                    }
                    break;
                case 2 :
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:29:22: ( '0' .. '9' )+ '.' ( '0' .. '9' )+
                    {
                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:29:22: ( '0' .. '9' )+
                    int cnt9=0;
                    loop9:
                    do {
                        int alt9=2;
                        int LA9_0 = input.LA(1);

                        if ( ((LA9_0 >= '0' && LA9_0 <= '9')) ) {
                            alt9=1;
                        }


                        switch (alt9) {
                    	case 1 :
                    	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
                    	        input.consume();
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;
                    	    }


                    	    }
                    	    break;

                    	default :
                    	    if ( cnt9 >= 1 ) break loop9;
                                EarlyExitException eee =
                                    new EarlyExitException(9, input);
                                throw eee;
                        }
                        cnt9++;
                    } while (true);


                    match('.'); 

                    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:29:36: ( '0' .. '9' )+
                    int cnt10=0;
                    loop10:
                    do {
                        int alt10=2;
                        int LA10_0 = input.LA(1);

                        if ( ((LA10_0 >= '0' && LA10_0 <= '9')) ) {
                            alt10=1;
                        }


                        switch (alt10) {
                    	case 1 :
                    	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
                    	        input.consume();
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        recover(mse);
                    	        throw mse;
                    	    }


                    	    }
                    	    break;

                    	default :
                    	    if ( cnt10 >= 1 ) break loop10;
                                EarlyExitException eee =
                                    new EarlyExitException(10, input);
                                throw eee;
                        }
                        cnt10++;
                    } while (true);


                    }
                    break;

            }
            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NUMBER"

    // $ANTLR start "COMMENT"
    public final void mCOMMENT() throws RecognitionException {
        try {
            int _type = COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:30:9: ( '/' '/' (~ NEWLINE )* ( NEWLINE )* )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:30:11: '/' '/' (~ NEWLINE )* ( NEWLINE )*
            {
            match('/'); 

            match('/'); 

            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:30:18: (~ NEWLINE )*
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( ((LA12_0 >= '\u0000' && LA12_0 <= '\t')||(LA12_0 >= '\u000B' && LA12_0 <= '\f')||(LA12_0 >= '\u000E' && LA12_0 <= '\uFFFF')) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\t')||(input.LA(1) >= '\u000B' && input.LA(1) <= '\f')||(input.LA(1) >= '\u000E' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop12;
                }
            } while (true);


            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:30:28: ( NEWLINE )*
            loop13:
            do {
                int alt13=2;
                int LA13_0 = input.LA(1);

                if ( (LA13_0=='\n'||LA13_0=='\r') ) {
                    alt13=1;
                }


                switch (alt13) {
            	case 1 :
            	    // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            	    {
            	    if ( input.LA(1)=='\n'||input.LA(1)=='\r' ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop13;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "COMMENT"

    // $ANTLR start "WHITESPACE"
    public final void mWHITESPACE() throws RecognitionException {
        try {
            int _type = WHITESPACE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:31:13: ( '\\t' | ' ' )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            {
            if ( input.LA(1)=='\t'||input.LA(1)==' ' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "WHITESPACE"

    // $ANTLR start "IDENTIFIER_CHAR"
    public final void mIDENTIFIER_CHAR() throws RecognitionException {
        try {
            int _type = IDENTIFIER_CHAR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:33:2: ( 'a' .. 'z' | 'A' .. 'Z' | '_' | '0' .. '9' | '%' | '-' )
            // /ews/QSPR-Core/src/qspr/oscript/OScript.g:
            {
            if ( input.LA(1)=='%'||input.LA(1)=='-'||(input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "IDENTIFIER_CHAR"

    public void mTokens() throws RecognitionException {
        // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:8: ( T__12 | T__13 | T__14 | IDENTIFIER | NEWLINE | ARROW | PREDICATE | NUMBER | COMMENT | WHITESPACE | IDENTIFIER_CHAR )
        int alt14=11;
        alt14 = dfa14.predict(input);
        switch (alt14) {
            case 1 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:10: T__12
                {
                mT__12(); 


                }
                break;
            case 2 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:16: T__13
                {
                mT__13(); 


                }
                break;
            case 3 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:22: T__14
                {
                mT__14(); 


                }
                break;
            case 4 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:28: IDENTIFIER
                {
                mIDENTIFIER(); 


                }
                break;
            case 5 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:39: NEWLINE
                {
                mNEWLINE(); 


                }
                break;
            case 6 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:47: ARROW
                {
                mARROW(); 


                }
                break;
            case 7 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:53: PREDICATE
                {
                mPREDICATE(); 


                }
                break;
            case 8 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:63: NUMBER
                {
                mNUMBER(); 


                }
                break;
            case 9 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:70: COMMENT
                {
                mCOMMENT(); 


                }
                break;
            case 10 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:78: WHITESPACE
                {
                mWHITESPACE(); 


                }
                break;
            case 11 :
                // /ews/QSPR-Core/src/qspr/oscript/OScript.g:1:89: IDENTIFIER_CHAR
                {
                mIDENTIFIER_CHAR(); 


                }
                break;

        }

    }


    protected DFA11 dfa11 = new DFA11(this);
    protected DFA14 dfa14 = new DFA14(this);
    static final String DFA11_eotS =
        "\1\uffff\1\2\2\uffff";
    static final String DFA11_eofS =
        "\4\uffff";
    static final String DFA11_minS =
        "\1\60\1\56\2\uffff";
    static final String DFA11_maxS =
        "\2\71\2\uffff";
    static final String DFA11_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA11_specialS =
        "\4\uffff}>";
    static final String[] DFA11_transitionS = {
            "\12\1",
            "\1\3\1\uffff\12\1",
            "",
            ""
    };

    static final short[] DFA11_eot = DFA.unpackEncodedString(DFA11_eotS);
    static final short[] DFA11_eof = DFA.unpackEncodedString(DFA11_eofS);
    static final char[] DFA11_min = DFA.unpackEncodedStringToUnsignedChars(DFA11_minS);
    static final char[] DFA11_max = DFA.unpackEncodedStringToUnsignedChars(DFA11_maxS);
    static final short[] DFA11_accept = DFA.unpackEncodedString(DFA11_acceptS);
    static final short[] DFA11_special = DFA.unpackEncodedString(DFA11_specialS);
    static final short[][] DFA11_transition;

    static {
        int numStates = DFA11_transitionS.length;
        DFA11_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA11_transition[i] = DFA.unpackEncodedString(DFA11_transitionS[i]);
        }
    }

    class DFA11 extends DFA {

        public DFA11(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 11;
            this.eot = DFA11_eot;
            this.eof = DFA11_eof;
            this.min = DFA11_min;
            this.max = DFA11_max;
            this.accept = DFA11_accept;
            this.special = DFA11_special;
            this.transition = DFA11_transition;
        }
        public String getDescription() {
            return "29:1: NUMBER : ( ( '0' .. '9' )+ | ( '0' .. '9' )+ '.' ( '0' .. '9' )+ );";
        }
    }
    static final String DFA14_eotS =
        "\7\uffff\1\17\1\14\10\uffff";
    static final String DFA14_eofS =
        "\21\uffff";
    static final String DFA14_minS =
        "\1\11\6\uffff\1\11\1\76\5\uffff\1\11\2\uffff";
    static final String DFA14_maxS =
        "\1\172\6\uffff\2\76\5\uffff\1\76\2\uffff";
    static final String DFA14_acceptS =
        "\1\uffff\1\1\1\2\1\3\2\4\1\5\2\uffff\1\7\1\10\1\11\1\13\1\6\1\uffff"+
        "\1\12\1\10";
    static final String DFA14_specialS =
        "\21\uffff}>";
    static final String[] DFA14_transitionS = {
            "\1\7\1\6\2\uffff\1\6\22\uffff\1\7\1\uffff\1\5\2\uffff\1\14\4"+
            "\uffff\1\1\2\uffff\1\10\1\uffff\1\13\12\12\2\uffff\3\11\2\uffff"+
            "\32\4\1\2\1\uffff\1\3\1\uffff\1\4\1\uffff\32\4",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\16\26\uffff\1\16\14\uffff\1\15\16\uffff\3\11",
            "\1\15",
            "",
            "",
            "",
            "",
            "",
            "\1\16\26\uffff\1\16\14\uffff\1\15\16\uffff\3\11",
            "",
            ""
    };

    static final short[] DFA14_eot = DFA.unpackEncodedString(DFA14_eotS);
    static final short[] DFA14_eof = DFA.unpackEncodedString(DFA14_eofS);
    static final char[] DFA14_min = DFA.unpackEncodedStringToUnsignedChars(DFA14_minS);
    static final char[] DFA14_max = DFA.unpackEncodedStringToUnsignedChars(DFA14_maxS);
    static final short[] DFA14_accept = DFA.unpackEncodedString(DFA14_acceptS);
    static final short[] DFA14_special = DFA.unpackEncodedString(DFA14_specialS);
    static final short[][] DFA14_transition;

    static {
        int numStates = DFA14_transitionS.length;
        DFA14_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA14_transition[i] = DFA.unpackEncodedString(DFA14_transitionS[i]);
        }
    }

    class DFA14 extends DFA {

        public DFA14(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 14;
            this.eot = DFA14_eot;
            this.eof = DFA14_eof;
            this.min = DFA14_min;
            this.max = DFA14_max;
            this.accept = DFA14_accept;
            this.special = DFA14_special;
            this.transition = DFA14_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__12 | T__13 | T__14 | IDENTIFIER | NEWLINE | ARROW | PREDICATE | NUMBER | COMMENT | WHITESPACE | IDENTIFIER_CHAR );";
        }
    }
 

}