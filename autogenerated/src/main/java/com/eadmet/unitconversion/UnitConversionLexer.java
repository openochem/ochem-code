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

package com.eadmet.unitconversion; 

import org.antlr.runtime.BaseRecognizer;
import org.antlr.runtime.CharStream;
import org.antlr.runtime.DFA;
import org.antlr.runtime.EarlyExitException;
import org.antlr.runtime.Lexer;
import org.antlr.runtime.MismatchedSetException;
import org.antlr.runtime.NoViableAltException;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.RecognizerSharedState;

public class UnitConversionLexer extends Lexer {
    public static final int E=8;
    public static final int LOG=19;
    public static final int NUMBER=25;
    public static final int WHITESPACE=26;
    public static final int FLOAT=28;
    public static final int INT=27;
    public static final int LOGN=21;
    public static final int MINUS=5;
    public static final int SQRT=13;
    public static final int EOF=-1;
    public static final int MUL=6;
    public static final int VARVALUE=23;
    public static final int LOGD=20;
    public static final int VARMW=24;
    public static final int LOGB=22;
    public static final int POW=14;
    public static final int LPAR=15;
    public static final int SIN=12;
    public static final int EXP=18;
    public static final int COMMA=16;
    public static final int COS=11;
    public static final int RPAR=17;
    public static final int PLUS=4;
    public static final int PI=9;
    public static final int DIGIT=29;
    public static final int DIV=7;
    public static final int NA=10;

    // delegates
    // delegators

    public UnitConversionLexer() {;} 
    public UnitConversionLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public UnitConversionLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);

    }
    public String getGrammarFileName() { return "UnitConversion.g"; }

    // $ANTLR start "WHITESPACE"
    public final void mWHITESPACE() throws RecognitionException {
        try {
            int _type = WHITESPACE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:92:12: ( ( '\\t' | ' ' | '\\r' | '\\n' | '\\u000C' )+ )
            // UnitConversion.g:92:14: ( '\\t' | ' ' | '\\r' | '\\n' | '\\u000C' )+
            {
            // UnitConversion.g:92:14: ( '\\t' | ' ' | '\\r' | '\\n' | '\\u000C' )+
            int cnt1=0;
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( ((LA1_0>='\t' && LA1_0<='\n')||(LA1_0>='\f' && LA1_0<='\r')||LA1_0==' ') ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // UnitConversion.g:
            	    {
            	    if ( (input.LA(1)>='\t' && input.LA(1)<='\n')||(input.LA(1)>='\f' && input.LA(1)<='\r')||input.LA(1)==' ' ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    if ( cnt1 >= 1 ) break loop1;
                        EarlyExitException eee =
                            new EarlyExitException(1, input);
                        throw eee;
                }
                cnt1++;
            } while (true);

             _channel = HIDDEN; 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "WHITESPACE"

    // $ANTLR start "NUMBER"
    public final void mNUMBER() throws RecognitionException {
        try {
            int _type = NUMBER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:96:9: ( INT | FLOAT )
            int alt2=2;
            alt2 = dfa2.predict(input);
            switch (alt2) {
                case 1 :
                    // UnitConversion.g:96:11: INT
                    {
                    mINT(); 

                    }
                    break;
                case 2 :
                    // UnitConversion.g:96:17: FLOAT
                    {
                    mFLOAT(); 

                    }
                    break;

            }
            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "NUMBER"

    // $ANTLR start "INT"
    public final void mINT() throws RecognitionException {
        try {
            // UnitConversion.g:98:14: ( ( DIGIT )+ )
            // UnitConversion.g:98:16: ( DIGIT )+
            {
            // UnitConversion.g:98:16: ( DIGIT )+
            int cnt3=0;
            loop3:
            do {
                int alt3=2;
                int LA3_0 = input.LA(1);

                if ( ((LA3_0>='0' && LA3_0<='9')) ) {
                    alt3=1;
                }


                switch (alt3) {
            	case 1 :
            	    // UnitConversion.g:98:16: DIGIT
            	    {
            	    mDIGIT(); 

            	    }
            	    break;

            	default :
            	    if ( cnt3 >= 1 ) break loop3;
                        EarlyExitException eee =
                            new EarlyExitException(3, input);
                        throw eee;
                }
                cnt3++;
            } while (true);


            }

        }
        finally {
        }
    }
    // $ANTLR end "INT"

    // $ANTLR start "FLOAT"
    public final void mFLOAT() throws RecognitionException {
        try {
            // UnitConversion.g:99:16: ( INT '.' ( DIGIT )* )
            // UnitConversion.g:99:18: INT '.' ( DIGIT )*
            {
            mINT(); 
            match('.'); 
            // UnitConversion.g:99:26: ( DIGIT )*
            loop4:
            do {
                int alt4=2;
                int LA4_0 = input.LA(1);

                if ( ((LA4_0>='0' && LA4_0<='9')) ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // UnitConversion.g:99:26: DIGIT
            	    {
            	    mDIGIT(); 

            	    }
            	    break;

            	default :
            	    break loop4;
                }
            } while (true);


            }

        }
        finally {
        }
    }
    // $ANTLR end "FLOAT"

    // $ANTLR start "DIGIT"
    public final void mDIGIT() throws RecognitionException {
        try {
            // UnitConversion.g:101:16: ( '0' .. '9' )
            // UnitConversion.g:101:18: '0' .. '9'
            {
            matchRange('0','9'); 

            }

        }
        finally {
        }
    }
    // $ANTLR end "DIGIT"

    // $ANTLR start "PLUS"
    public final void mPLUS() throws RecognitionException {
        try {
            int _type = PLUS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:105:9: ( '+' )
            // UnitConversion.g:105:11: '+'
            {
            match('+'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "PLUS"

    // $ANTLR start "MINUS"
    public final void mMINUS() throws RecognitionException {
        try {
            int _type = MINUS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:106:9: ( '-' )
            // UnitConversion.g:106:11: '-'
            {
            match('-'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "MINUS"

    // $ANTLR start "MUL"
    public final void mMUL() throws RecognitionException {
        try {
            int _type = MUL;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:107:9: ( '*' )
            // UnitConversion.g:107:11: '*'
            {
            match('*'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "MUL"

    // $ANTLR start "DIV"
    public final void mDIV() throws RecognitionException {
        try {
            int _type = DIV;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:108:9: ( '/' )
            // UnitConversion.g:108:11: '/'
            {
            match('/'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "DIV"

    // $ANTLR start "COMMA"
    public final void mCOMMA() throws RecognitionException {
        try {
            int _type = COMMA;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:112:9: ( ',' )
            // UnitConversion.g:112:11: ','
            {
            match(','); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "COMMA"

    // $ANTLR start "COS"
    public final void mCOS() throws RecognitionException {
        try {
            int _type = COS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:114:9: ( 'cos' )
            // UnitConversion.g:114:11: 'cos'
            {
            match("cos"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "COS"

    // $ANTLR start "SIN"
    public final void mSIN() throws RecognitionException {
        try {
            int _type = SIN;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:115:9: ( 'sin' )
            // UnitConversion.g:115:11: 'sin'
            {
            match("sin"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "SIN"

    // $ANTLR start "SQRT"
    public final void mSQRT() throws RecognitionException {
        try {
            int _type = SQRT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:116:9: ( 'sqrt' )
            // UnitConversion.g:116:11: 'sqrt'
            {
            match("sqrt"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "SQRT"

    // $ANTLR start "POW"
    public final void mPOW() throws RecognitionException {
        try {
            int _type = POW;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:117:9: ( 'pow' )
            // UnitConversion.g:117:11: 'pow'
            {
            match("pow"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "POW"

    // $ANTLR start "EXP"
    public final void mEXP() throws RecognitionException {
        try {
            int _type = EXP;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:118:9: ( 'exp' )
            // UnitConversion.g:118:11: 'exp'
            {
            match("exp"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "EXP"

    // $ANTLR start "LOG"
    public final void mLOG() throws RecognitionException {
        try {
            int _type = LOG;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:119:9: ( 'log' )
            // UnitConversion.g:119:11: 'log'
            {
            match("log"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "LOG"

    // $ANTLR start "LOGD"
    public final void mLOGD() throws RecognitionException {
        try {
            int _type = LOGD;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:120:9: ( 'lg' | 'log10' )
            int alt5=2;
            int LA5_0 = input.LA(1);

            if ( (LA5_0=='l') ) {
                int LA5_1 = input.LA(2);

                if ( (LA5_1=='g') ) {
                    alt5=1;
                }
                else if ( (LA5_1=='o') ) {
                    alt5=2;
                }
                else {
                    NoViableAltException nvae =
                        new NoViableAltException("", 5, 1, input);

                    throw nvae;
                }
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;
            }
            switch (alt5) {
                case 1 :
                    // UnitConversion.g:120:11: 'lg'
                    {
                    match("lg"); 


                    }
                    break;
                case 2 :
                    // UnitConversion.g:120:18: 'log10'
                    {
                    match("log10"); 


                    }
                    break;

            }
            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "LOGD"

    // $ANTLR start "LOGN"
    public final void mLOGN() throws RecognitionException {
        try {
            int _type = LOGN;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:121:9: ( 'ln' | 'logn' )
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0=='l') ) {
                int LA6_1 = input.LA(2);

                if ( (LA6_1=='n') ) {
                    alt6=1;
                }
                else if ( (LA6_1=='o') ) {
                    alt6=2;
                }
                else {
                    NoViableAltException nvae =
                        new NoViableAltException("", 6, 1, input);

                    throw nvae;
                }
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 6, 0, input);

                throw nvae;
            }
            switch (alt6) {
                case 1 :
                    // UnitConversion.g:121:11: 'ln'
                    {
                    match("ln"); 


                    }
                    break;
                case 2 :
                    // UnitConversion.g:121:18: 'logn'
                    {
                    match("logn"); 


                    }
                    break;

            }
            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "LOGN"

    // $ANTLR start "LOGB"
    public final void mLOGB() throws RecognitionException {
        try {
            int _type = LOGB;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:122:9: ( 'log2' )
            // UnitConversion.g:122:11: 'log2'
            {
            match("log2"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "LOGB"

    // $ANTLR start "E"
    public final void mE() throws RecognitionException {
        try {
            int _type = E;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:126:9: ( 'E' )
            // UnitConversion.g:126:11: 'E'
            {
            match('E'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "E"

    // $ANTLR start "PI"
    public final void mPI() throws RecognitionException {
        try {
            int _type = PI;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:127:9: ( 'PI' )
            // UnitConversion.g:127:11: 'PI'
            {
            match("PI"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "PI"

    // $ANTLR start "NA"
    public final void mNA() throws RecognitionException {
        try {
            int _type = NA;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:128:9: ( 'NA' )
            // UnitConversion.g:128:11: 'NA'
            {
            match("NA"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "NA"

    // $ANTLR start "LPAR"
    public final void mLPAR() throws RecognitionException {
        try {
            int _type = LPAR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:132:6: ( '(' )
            // UnitConversion.g:132:8: '('
            {
            match('('); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "LPAR"

    // $ANTLR start "RPAR"
    public final void mRPAR() throws RecognitionException {
        try {
            int _type = RPAR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:133:6: ( ')' )
            // UnitConversion.g:133:8: ')'
            {
            match(')'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "RPAR"

    // $ANTLR start "VARVALUE"
    public final void mVARVALUE() throws RecognitionException {
        try {
            int _type = VARVALUE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:137:10: ( '#' )
            // UnitConversion.g:137:12: '#'
            {
            match('#'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "VARVALUE"

    // $ANTLR start "VARMW"
    public final void mVARMW() throws RecognitionException {
        try {
            int _type = VARMW;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // UnitConversion.g:138:10: ( '$MW' )
            // UnitConversion.g:138:12: '$MW'
            {
            match("$MW"); 


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "VARMW"

    public void mTokens() throws RecognitionException {
        // UnitConversion.g:1:8: ( WHITESPACE | NUMBER | PLUS | MINUS | MUL | DIV | COMMA | COS | SIN | SQRT | POW | EXP | LOG | LOGD | LOGN | LOGB | E | PI | NA | LPAR | RPAR | VARVALUE | VARMW )
        int alt7=23;
        alt7 = dfa7.predict(input);
        switch (alt7) {
            case 1 :
                // UnitConversion.g:1:10: WHITESPACE
                {
                mWHITESPACE(); 

                }
                break;
            case 2 :
                // UnitConversion.g:1:21: NUMBER
                {
                mNUMBER(); 

                }
                break;
            case 3 :
                // UnitConversion.g:1:28: PLUS
                {
                mPLUS(); 

                }
                break;
            case 4 :
                // UnitConversion.g:1:33: MINUS
                {
                mMINUS(); 

                }
                break;
            case 5 :
                // UnitConversion.g:1:39: MUL
                {
                mMUL(); 

                }
                break;
            case 6 :
                // UnitConversion.g:1:43: DIV
                {
                mDIV(); 

                }
                break;
            case 7 :
                // UnitConversion.g:1:47: COMMA
                {
                mCOMMA(); 

                }
                break;
            case 8 :
                // UnitConversion.g:1:53: COS
                {
                mCOS(); 

                }
                break;
            case 9 :
                // UnitConversion.g:1:57: SIN
                {
                mSIN(); 

                }
                break;
            case 10 :
                // UnitConversion.g:1:61: SQRT
                {
                mSQRT(); 

                }
                break;
            case 11 :
                // UnitConversion.g:1:66: POW
                {
                mPOW(); 

                }
                break;
            case 12 :
                // UnitConversion.g:1:70: EXP
                {
                mEXP(); 

                }
                break;
            case 13 :
                // UnitConversion.g:1:74: LOG
                {
                mLOG(); 

                }
                break;
            case 14 :
                // UnitConversion.g:1:78: LOGD
                {
                mLOGD(); 

                }
                break;
            case 15 :
                // UnitConversion.g:1:83: LOGN
                {
                mLOGN(); 

                }
                break;
            case 16 :
                // UnitConversion.g:1:88: LOGB
                {
                mLOGB(); 

                }
                break;
            case 17 :
                // UnitConversion.g:1:93: E
                {
                mE(); 

                }
                break;
            case 18 :
                // UnitConversion.g:1:95: PI
                {
                mPI(); 

                }
                break;
            case 19 :
                // UnitConversion.g:1:98: NA
                {
                mNA(); 

                }
                break;
            case 20 :
                // UnitConversion.g:1:101: LPAR
                {
                mLPAR(); 

                }
                break;
            case 21 :
                // UnitConversion.g:1:106: RPAR
                {
                mRPAR(); 

                }
                break;
            case 22 :
                // UnitConversion.g:1:111: VARVALUE
                {
                mVARVALUE(); 

                }
                break;
            case 23 :
                // UnitConversion.g:1:120: VARMW
                {
                mVARMW(); 

                }
                break;

        }

    }


    protected DFA2 dfa2 = new DFA2(this);
    protected DFA7 dfa7 = new DFA7(this);
    static final String DFA2_eotS =
        "\1\uffff\1\2\2\uffff";
    static final String DFA2_eofS =
        "\4\uffff";
    static final String DFA2_minS =
        "\1\60\1\56\2\uffff";
    static final String DFA2_maxS =
        "\2\71\2\uffff";
    static final String DFA2_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA2_specialS =
        "\4\uffff}>";
    static final String[] DFA2_transitionS = {
            "\12\1",
            "\1\3\1\uffff\12\1",
            "",
            ""
    };

    static final short[] DFA2_eot = DFA.unpackEncodedString(DFA2_eotS);
    static final short[] DFA2_eof = DFA.unpackEncodedString(DFA2_eofS);
    static final char[] DFA2_min = DFA.unpackEncodedStringToUnsignedChars(DFA2_minS);
    static final char[] DFA2_max = DFA.unpackEncodedStringToUnsignedChars(DFA2_maxS);
    static final short[] DFA2_accept = DFA.unpackEncodedString(DFA2_acceptS);
    static final short[] DFA2_special = DFA.unpackEncodedString(DFA2_specialS);
    static final short[][] DFA2_transition;

    static {
        int numStates = DFA2_transitionS.length;
        DFA2_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA2_transition[i] = DFA.unpackEncodedString(DFA2_transitionS[i]);
        }
    }

    class DFA2 extends DFA {

        public DFA2(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 2;
            this.eot = DFA2_eot;
            this.eof = DFA2_eof;
            this.min = DFA2_min;
            this.max = DFA2_max;
            this.accept = DFA2_accept;
            this.special = DFA2_special;
            this.transition = DFA2_transition;
        }
        public String getDescription() {
            return "96:1: NUMBER : ( INT | FLOAT );";
        }
    }
    static final String DFA7_eotS =
        "\31\uffff\1\33\2\uffff";
    static final String DFA7_eofS =
        "\34\uffff";
    static final String DFA7_minS =
        "\1\11\10\uffff\1\151\2\uffff\1\147\11\uffff\1\147\2\uffff\1\61\2"+
        "\uffff";
    static final String DFA7_maxS =
        "\1\163\10\uffff\1\161\2\uffff\1\157\11\uffff\1\147\2\uffff\1\156"+
        "\2\uffff";
    static final String DFA7_acceptS =
        "\1\uffff\1\1\1\2\1\3\1\4\1\5\1\6\1\7\1\10\1\uffff\1\13\1\14\1\uffff"+
        "\1\21\1\22\1\23\1\24\1\25\1\26\1\27\1\11\1\12\1\uffff\1\16\1\17"+
        "\1\uffff\1\20\1\15";
    static final String DFA7_specialS =
        "\34\uffff}>";
    static final String[] DFA7_transitionS = {
            "\2\1\1\uffff\2\1\22\uffff\1\1\2\uffff\1\22\1\23\3\uffff\1\20"+
            "\1\21\1\5\1\3\1\7\1\4\1\uffff\1\6\12\2\13\uffff\1\15\10\uffff"+
            "\1\17\1\uffff\1\16\22\uffff\1\10\1\uffff\1\13\6\uffff\1\14\3"+
            "\uffff\1\12\2\uffff\1\11",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\24\7\uffff\1\25",
            "",
            "",
            "\1\27\6\uffff\1\30\1\26",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\31",
            "",
            "",
            "\1\27\1\32\73\uffff\1\30",
            "",
            ""
    };

    static final short[] DFA7_eot = DFA.unpackEncodedString(DFA7_eotS);
    static final short[] DFA7_eof = DFA.unpackEncodedString(DFA7_eofS);
    static final char[] DFA7_min = DFA.unpackEncodedStringToUnsignedChars(DFA7_minS);
    static final char[] DFA7_max = DFA.unpackEncodedStringToUnsignedChars(DFA7_maxS);
    static final short[] DFA7_accept = DFA.unpackEncodedString(DFA7_acceptS);
    static final short[] DFA7_special = DFA.unpackEncodedString(DFA7_specialS);
    static final short[][] DFA7_transition;

    static {
        int numStates = DFA7_transitionS.length;
        DFA7_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA7_transition[i] = DFA.unpackEncodedString(DFA7_transitionS[i]);
        }
    }

    class DFA7 extends DFA {

        public DFA7(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 7;
            this.eot = DFA7_eot;
            this.eof = DFA7_eof;
            this.min = DFA7_min;
            this.max = DFA7_max;
            this.accept = DFA7_accept;
            this.special = DFA7_special;
            this.transition = DFA7_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( WHITESPACE | NUMBER | PLUS | MINUS | MUL | DIV | COMMA | COS | SIN | SQRT | POW | EXP | LOG | LOGD | LOGN | LOGB | E | PI | NA | LPAR | RPAR | VARVALUE | VARMW );";
        }
    }
 

}