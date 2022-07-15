// $ANTLR 3.1.3 Mar 18, 2009 10:09:25 UnitConversion.g 2009-05-27 14:26:42
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

import org.antlr.runtime.BitSet;
import org.antlr.runtime.NoViableAltException;
import org.antlr.runtime.Parser;
import org.antlr.runtime.RecognitionException;
import org.antlr.runtime.RecognizerSharedState;
import org.antlr.runtime.Token;
import org.antlr.runtime.TokenStream;

public class UnitConversionParser extends Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "PLUS", "MINUS", "MUL", "DIV", "E", "PI", "NA", "COS", "SIN", "SQRT", "POW", "LPAR", "COMMA", "RPAR", "EXP", "LOG", "LOGD", "LOGN", "LOGB", "VARVALUE", "VARMW", "NUMBER", "WHITESPACE", "INT", "FLOAT", "DIGIT"
    };
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
    public static final int PLUS=4;
    public static final int RPAR=17;
    public static final int PI=9;
    public static final int DIGIT=29;
    public static final int DIV=7;
    public static final int NA=10;

    // delegates
    // delegators


        public UnitConversionParser(TokenStream input) {
            this(input, new RecognizerSharedState());
        }
        public UnitConversionParser(TokenStream input, RecognizerSharedState state) {
            super(input, state);
             
        }
        

    public String[] getTokenNames() { return UnitConversionParser.tokenNames; }
    public String getGrammarFileName() { return "UnitConversion.g"; }


        public double varvalue;  // Value to-be-converted.
        public double varmw;     // Molecular weight of compound. 
        
        public static double log(double base, double arg) { return Math.log(arg) / Math.log(base); }  // Logarithm with arbitrary base.



    // $ANTLR start "start"
    // UnitConversion.g:19:1: start returns [double value] : sumexpr EOF ;
    public final double start() throws RecognitionException {
        double value = 0.0;

        double sumexpr1 = 0.0;


        try {
            // UnitConversion.g:20:5: ( sumexpr EOF )
            // UnitConversion.g:20:7: sumexpr EOF
            {
            pushFollow(FOLLOW_sumexpr_in_start43);
            sumexpr1=sumexpr();

            state._fsp--;

             value = sumexpr1; 
            match(input,EOF,FOLLOW_EOF_in_start47); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "start"


    // $ANTLR start "sumexpr"
    // UnitConversion.g:23:1: sumexpr returns [double value] : v= mulexpr ( ( PLUS w= mulexpr ) | ( MINUS w= mulexpr ) )* ;
    public final double sumexpr() throws RecognitionException {
        double value = 0.0;

        double v = 0.0;

        double w = 0.0;


        try {
            // UnitConversion.g:24:5: (v= mulexpr ( ( PLUS w= mulexpr ) | ( MINUS w= mulexpr ) )* )
            // UnitConversion.g:24:7: v= mulexpr ( ( PLUS w= mulexpr ) | ( MINUS w= mulexpr ) )*
            {
            pushFollow(FOLLOW_mulexpr_in_sumexpr67);
            v=mulexpr();

            state._fsp--;

             value = v; 
            // UnitConversion.g:25:11: ( ( PLUS w= mulexpr ) | ( MINUS w= mulexpr ) )*
            loop1:
            do {
                int alt1=3;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==PLUS) ) {
                    alt1=1;
                }
                else if ( (LA1_0==MINUS) ) {
                    alt1=2;
                }


                switch (alt1) {
            	case 1 :
            	    // UnitConversion.g:26:13: ( PLUS w= mulexpr )
            	    {
            	    // UnitConversion.g:26:13: ( PLUS w= mulexpr )
            	    // UnitConversion.g:26:15: PLUS w= mulexpr
            	    {
            	    match(input,PLUS,FOLLOW_PLUS_in_sumexpr98); 
            	    pushFollow(FOLLOW_mulexpr_in_sumexpr103);
            	    w=mulexpr();

            	    state._fsp--;

            	     value += w; 

            	    }


            	    }
            	    break;
            	case 2 :
            	    // UnitConversion.g:27:13: ( MINUS w= mulexpr )
            	    {
            	    // UnitConversion.g:27:13: ( MINUS w= mulexpr )
            	    // UnitConversion.g:27:15: MINUS w= mulexpr
            	    {
            	    match(input,MINUS,FOLLOW_MINUS_in_sumexpr124); 
            	    pushFollow(FOLLOW_mulexpr_in_sumexpr128);
            	    w=mulexpr();

            	    state._fsp--;

            	     value -= w; 

            	    }


            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "sumexpr"


    // $ANTLR start "mulexpr"
    // UnitConversion.g:31:1: mulexpr returns [double value] : v= unaryexpr ( ( MUL w= unaryexpr ) | ( DIV w= unaryexpr ) )* ;
    public final double mulexpr() throws RecognitionException {
        double value = 0.0;

        double v = 0.0;

        double w = 0.0;


        try {
            // UnitConversion.g:32:5: (v= unaryexpr ( ( MUL w= unaryexpr ) | ( DIV w= unaryexpr ) )* )
            // UnitConversion.g:32:7: v= unaryexpr ( ( MUL w= unaryexpr ) | ( DIV w= unaryexpr ) )*
            {
            pushFollow(FOLLOW_unaryexpr_in_mulexpr166);
            v=unaryexpr();

            state._fsp--;

             value = v; 
            // UnitConversion.g:33:11: ( ( MUL w= unaryexpr ) | ( DIV w= unaryexpr ) )*
            loop2:
            do {
                int alt2=3;
                int LA2_0 = input.LA(1);

                if ( (LA2_0==MUL) ) {
                    alt2=1;
                }
                else if ( (LA2_0==DIV) ) {
                    alt2=2;
                }


                switch (alt2) {
            	case 1 :
            	    // UnitConversion.g:34:13: ( MUL w= unaryexpr )
            	    {
            	    // UnitConversion.g:34:13: ( MUL w= unaryexpr )
            	    // UnitConversion.g:34:15: MUL w= unaryexpr
            	    {
            	    match(input,MUL,FOLLOW_MUL_in_mulexpr197); 
            	    pushFollow(FOLLOW_unaryexpr_in_mulexpr202);
            	    w=unaryexpr();

            	    state._fsp--;

            	     value *= w; 

            	    }


            	    }
            	    break;
            	case 2 :
            	    // UnitConversion.g:35:13: ( DIV w= unaryexpr )
            	    {
            	    // UnitConversion.g:35:13: ( DIV w= unaryexpr )
            	    // UnitConversion.g:35:15: DIV w= unaryexpr
            	    {
            	    match(input,DIV,FOLLOW_DIV_in_mulexpr223); 
            	    pushFollow(FOLLOW_unaryexpr_in_mulexpr228);
            	    w=unaryexpr();

            	    state._fsp--;

            	     value /= w; 

            	    }


            	    }
            	    break;

            	default :
            	    break loop2;
                }
            } while (true);


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "mulexpr"


    // $ANTLR start "unaryexpr"
    // UnitConversion.g:39:1: unaryexpr returns [double value] : ( ( MINUS primexpr ) | ( primexpr ) );
    public final double unaryexpr() throws RecognitionException {
        double value = 0.0;

        double primexpr2 = 0.0;

        double primexpr3 = 0.0;


        try {
            // UnitConversion.g:40:5: ( ( MINUS primexpr ) | ( primexpr ) )
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( (LA3_0==MINUS) ) {
                alt3=1;
            }
            else if ( ((LA3_0>=E && LA3_0<=LPAR)||(LA3_0>=EXP && LA3_0<=NUMBER)) ) {
                alt3=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 3, 0, input);

                throw nvae;
            }
            switch (alt3) {
                case 1 :
                    // UnitConversion.g:41:11: ( MINUS primexpr )
                    {
                    // UnitConversion.g:41:11: ( MINUS primexpr )
                    // UnitConversion.g:41:13: MINUS primexpr
                    {
                    match(input,MINUS,FOLLOW_MINUS_in_unaryexpr279); 
                    pushFollow(FOLLOW_primexpr_in_unaryexpr281);
                    primexpr2=primexpr();

                    state._fsp--;

                     value = -primexpr2; 

                    }


                    }
                    break;
                case 2 :
                    // UnitConversion.g:42:11: ( primexpr )
                    {
                    // UnitConversion.g:42:11: ( primexpr )
                    // UnitConversion.g:42:13: primexpr
                    {
                    pushFollow(FOLLOW_primexpr_in_unaryexpr299);
                    primexpr3=primexpr();

                    state._fsp--;

                     value = primexpr3; 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "unaryexpr"


    // $ANTLR start "primexpr"
    // UnitConversion.g:46:1: primexpr returns [double value] : ( ( number ) | ( var ) | ( constant ) | ( funexpr ) | ( parexpr ) );
    public final double primexpr() throws RecognitionException {
        double value = 0.0;

        double number4 = 0.0;

        double var5 = 0.0;

        double constant6 = 0.0;

        double funexpr7 = 0.0;

        double parexpr8 = 0.0;


        try {
            // UnitConversion.g:47:5: ( ( number ) | ( var ) | ( constant ) | ( funexpr ) | ( parexpr ) )
            int alt4=5;
            switch ( input.LA(1) ) {
            case NUMBER:
                {
                alt4=1;
                }
                break;
            case VARVALUE:
            case VARMW:
                {
                alt4=2;
                }
                break;
            case E:
            case PI:
            case NA:
                {
                alt4=3;
                }
                break;
            case COS:
            case SIN:
            case SQRT:
            case POW:
            case EXP:
            case LOG:
            case LOGD:
            case LOGN:
            case LOGB:
                {
                alt4=4;
                }
                break;
            case LPAR:
                {
                alt4=5;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 4, 0, input);

                throw nvae;
            }

            switch (alt4) {
                case 1 :
                    // UnitConversion.g:47:11: ( number )
                    {
                    // UnitConversion.g:47:11: ( number )
                    // UnitConversion.g:47:13: number
                    {
                    pushFollow(FOLLOW_number_in_primexpr335);
                    number4=number();

                    state._fsp--;

                     value = number4;  

                    }


                    }
                    break;
                case 2 :
                    // UnitConversion.g:48:11: ( var )
                    {
                    // UnitConversion.g:48:11: ( var )
                    // UnitConversion.g:48:13: var
                    {
                    pushFollow(FOLLOW_var_in_primexpr355);
                    var5=var();

                    state._fsp--;

                     value = var5;     

                    }


                    }
                    break;
                case 3 :
                    // UnitConversion.g:49:11: ( constant )
                    {
                    // UnitConversion.g:49:11: ( constant )
                    // UnitConversion.g:49:13: constant
                    {
                    pushFollow(FOLLOW_constant_in_primexpr378);
                    constant6=constant();

                    state._fsp--;

                     value = constant6;

                    }


                    }
                    break;
                case 4 :
                    // UnitConversion.g:50:11: ( funexpr )
                    {
                    // UnitConversion.g:50:11: ( funexpr )
                    // UnitConversion.g:50:13: funexpr
                    {
                    pushFollow(FOLLOW_funexpr_in_primexpr396);
                    funexpr7=funexpr();

                    state._fsp--;

                     value = funexpr7; 

                    }


                    }
                    break;
                case 5 :
                    // UnitConversion.g:51:11: ( parexpr )
                    {
                    // UnitConversion.g:51:11: ( parexpr )
                    // UnitConversion.g:51:13: parexpr
                    {
                    pushFollow(FOLLOW_parexpr_in_primexpr415);
                    parexpr8=parexpr();

                    state._fsp--;

                     value = parexpr8; 

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "primexpr"


    // $ANTLR start "constant"
    // UnitConversion.g:55:1: constant returns [double value] : ( ( E ) | ( PI ) | ( NA ) );
    public final double constant() throws RecognitionException {
        double value = 0.0;

        try {
            // UnitConversion.g:56:5: ( ( E ) | ( PI ) | ( NA ) )
            int alt5=3;
            switch ( input.LA(1) ) {
            case E:
                {
                alt5=1;
                }
                break;
            case PI:
                {
                alt5=2;
                }
                break;
            case NA:
                {
                alt5=3;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;
            }

            switch (alt5) {
                case 1 :
                    // UnitConversion.g:56:11: ( E )
                    {
                    // UnitConversion.g:56:11: ( E )
                    // UnitConversion.g:56:13: E
                    {
                    match(input,E,FOLLOW_E_in_constant455); 
                     value = 2.7182818284590452; 

                    }


                    }
                    break;
                case 2 :
                    // UnitConversion.g:57:10: ( PI )
                    {
                    // UnitConversion.g:57:10: ( PI )
                    // UnitConversion.g:57:12: PI
                    {
                    match(input,PI,FOLLOW_PI_in_constant473); 
                     value = 3.1415926535897932; 

                    }


                    }
                    break;
                case 3 :
                    // UnitConversion.g:58:10: ( NA )
                    {
                    // UnitConversion.g:58:10: ( NA )
                    // UnitConversion.g:58:12: NA
                    {
                    match(input,NA,FOLLOW_NA_in_constant490); 
                     value = 6.02214179e23;      

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "constant"


    // $ANTLR start "funexpr"
    // UnitConversion.g:61:1: funexpr returns [double value] : ( ( COS parexpr ) | ( SIN parexpr ) | ( SQRT parexpr ) | ( POW LPAR v= sumexpr COMMA w= sumexpr RPAR ) | ( EXP parexpr ) | ( LOG LPAR v= sumexpr COMMA w= sumexpr RPAR ) | ( LOGD parexpr ) | ( LOGN parexpr ) | ( LOGB parexpr ) );
    public final double funexpr() throws RecognitionException {
        double value = 0.0;

        double v = 0.0;

        double w = 0.0;

        double parexpr9 = 0.0;

        double parexpr10 = 0.0;

        double parexpr11 = 0.0;

        double parexpr12 = 0.0;

        double parexpr13 = 0.0;

        double parexpr14 = 0.0;

        double parexpr15 = 0.0;


        try {
            // UnitConversion.g:62:5: ( ( COS parexpr ) | ( SIN parexpr ) | ( SQRT parexpr ) | ( POW LPAR v= sumexpr COMMA w= sumexpr RPAR ) | ( EXP parexpr ) | ( LOG LPAR v= sumexpr COMMA w= sumexpr RPAR ) | ( LOGD parexpr ) | ( LOGN parexpr ) | ( LOGB parexpr ) )
            int alt6=9;
            switch ( input.LA(1) ) {
            case COS:
                {
                alt6=1;
                }
                break;
            case SIN:
                {
                alt6=2;
                }
                break;
            case SQRT:
                {
                alt6=3;
                }
                break;
            case POW:
                {
                alt6=4;
                }
                break;
            case EXP:
                {
                alt6=5;
                }
                break;
            case LOG:
                {
                alt6=6;
                }
                break;
            case LOGD:
                {
                alt6=7;
                }
                break;
            case LOGN:
                {
                alt6=8;
                }
                break;
            case LOGB:
                {
                alt6=9;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 6, 0, input);

                throw nvae;
            }

            switch (alt6) {
                case 1 :
                    // UnitConversion.g:62:11: ( COS parexpr )
                    {
                    // UnitConversion.g:62:11: ( COS parexpr )
                    // UnitConversion.g:62:13: COS parexpr
                    {
                    match(input,COS,FOLLOW_COS_in_funexpr524); 
                    pushFollow(FOLLOW_parexpr_in_funexpr527);
                    parexpr9=parexpr();

                    state._fsp--;

                     value = Math.cos  (parexpr9); 

                    }


                    }
                    break;
                case 2 :
                    // UnitConversion.g:63:10: ( SIN parexpr )
                    {
                    // UnitConversion.g:63:10: ( SIN parexpr )
                    // UnitConversion.g:63:12: SIN parexpr
                    {
                    match(input,SIN,FOLLOW_SIN_in_funexpr544); 
                    pushFollow(FOLLOW_parexpr_in_funexpr547);
                    parexpr10=parexpr();

                    state._fsp--;

                     value = Math.sin  (parexpr10); 

                    }


                    }
                    break;
                case 3 :
                    // UnitConversion.g:64:10: ( SQRT parexpr )
                    {
                    // UnitConversion.g:64:10: ( SQRT parexpr )
                    // UnitConversion.g:64:12: SQRT parexpr
                    {
                    match(input,SQRT,FOLLOW_SQRT_in_funexpr564); 
                    pushFollow(FOLLOW_parexpr_in_funexpr566);
                    parexpr11=parexpr();

                    state._fsp--;

                     value = Math.sqrt (parexpr11); 

                    }


                    }
                    break;
                case 4 :
                    // UnitConversion.g:65:10: ( POW LPAR v= sumexpr COMMA w= sumexpr RPAR )
                    {
                    // UnitConversion.g:65:10: ( POW LPAR v= sumexpr COMMA w= sumexpr RPAR )
                    // UnitConversion.g:65:12: POW LPAR v= sumexpr COMMA w= sumexpr RPAR
                    {
                    match(input,POW,FOLLOW_POW_in_funexpr583); 
                    match(input,LPAR,FOLLOW_LPAR_in_funexpr586); 
                    pushFollow(FOLLOW_sumexpr_in_funexpr590);
                    v=sumexpr();

                    state._fsp--;

                    match(input,COMMA,FOLLOW_COMMA_in_funexpr592); 
                    pushFollow(FOLLOW_sumexpr_in_funexpr596);
                    w=sumexpr();

                    state._fsp--;

                    match(input,RPAR,FOLLOW_RPAR_in_funexpr598); 
                     value = Math.pow(v,w); 

                    }


                    }
                    break;
                case 5 :
                    // UnitConversion.g:66:10: ( EXP parexpr )
                    {
                    // UnitConversion.g:66:10: ( EXP parexpr )
                    // UnitConversion.g:66:12: EXP parexpr
                    {
                    match(input,EXP,FOLLOW_EXP_in_funexpr615); 
                    pushFollow(FOLLOW_parexpr_in_funexpr618);
                    parexpr12=parexpr();

                    state._fsp--;

                     value = Math.exp  (parexpr12); 

                    }


                    }
                    break;
                case 6 :
                    // UnitConversion.g:67:10: ( LOG LPAR v= sumexpr COMMA w= sumexpr RPAR )
                    {
                    // UnitConversion.g:67:10: ( LOG LPAR v= sumexpr COMMA w= sumexpr RPAR )
                    // UnitConversion.g:67:12: LOG LPAR v= sumexpr COMMA w= sumexpr RPAR
                    {
                    match(input,LOG,FOLLOW_LOG_in_funexpr635); 
                    match(input,LPAR,FOLLOW_LPAR_in_funexpr638); 
                    pushFollow(FOLLOW_sumexpr_in_funexpr642);
                    v=sumexpr();

                    state._fsp--;

                    match(input,COMMA,FOLLOW_COMMA_in_funexpr644); 
                    pushFollow(FOLLOW_sumexpr_in_funexpr648);
                    w=sumexpr();

                    state._fsp--;

                    match(input,RPAR,FOLLOW_RPAR_in_funexpr650); 
                     value = log(v,w); 

                    }


                    }
                    break;
                case 7 :
                    // UnitConversion.g:68:10: ( LOGD parexpr )
                    {
                    // UnitConversion.g:68:10: ( LOGD parexpr )
                    // UnitConversion.g:68:12: LOGD parexpr
                    {
                    match(input,LOGD,FOLLOW_LOGD_in_funexpr667); 
                    pushFollow(FOLLOW_parexpr_in_funexpr669);
                    parexpr13=parexpr();

                    state._fsp--;

                     value = Math.log10(parexpr13); 

                    }


                    }
                    break;
                case 8 :
                    // UnitConversion.g:69:10: ( LOGN parexpr )
                    {
                    // UnitConversion.g:69:10: ( LOGN parexpr )
                    // UnitConversion.g:69:12: LOGN parexpr
                    {
                    match(input,LOGN,FOLLOW_LOGN_in_funexpr686); 
                    pushFollow(FOLLOW_parexpr_in_funexpr688);
                    parexpr14=parexpr();

                    state._fsp--;

                     value = Math.log  (parexpr14); 

                    }


                    }
                    break;
                case 9 :
                    // UnitConversion.g:70:10: ( LOGB parexpr )
                    {
                    // UnitConversion.g:70:10: ( LOGB parexpr )
                    // UnitConversion.g:70:12: LOGB parexpr
                    {
                    match(input,LOGB,FOLLOW_LOGB_in_funexpr705); 
                    pushFollow(FOLLOW_parexpr_in_funexpr707);
                    parexpr15=parexpr();

                    state._fsp--;

                     value = log(2.,parexpr15);     

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "funexpr"


    // $ANTLR start "parexpr"
    // UnitConversion.g:73:1: parexpr returns [double value] : LPAR sumexpr RPAR ;
    public final double parexpr() throws RecognitionException {
        double value = 0.0;

        double sumexpr16 = 0.0;


        try {
            // UnitConversion.g:74:5: ( LPAR sumexpr RPAR )
            // UnitConversion.g:74:7: LPAR sumexpr RPAR
            {
            match(input,LPAR,FOLLOW_LPAR_in_parexpr736); 
            pushFollow(FOLLOW_sumexpr_in_parexpr738);
            sumexpr16=sumexpr();

            state._fsp--;

             value = sumexpr16; 
            match(input,RPAR,FOLLOW_RPAR_in_parexpr742); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "parexpr"


    // $ANTLR start "var"
    // UnitConversion.g:77:1: var returns [double value] : ( ( VARVALUE ) | ( VARMW ) );
    public final double var() throws RecognitionException {
        double value = 0.0;

        try {
            // UnitConversion.g:78:5: ( ( VARVALUE ) | ( VARMW ) )
            int alt7=2;
            int LA7_0 = input.LA(1);

            if ( (LA7_0==VARVALUE) ) {
                alt7=1;
            }
            else if ( (LA7_0==VARMW) ) {
                alt7=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 7, 0, input);

                throw nvae;
            }
            switch (alt7) {
                case 1 :
                    // UnitConversion.g:78:10: ( VARVALUE )
                    {
                    // UnitConversion.g:78:10: ( VARVALUE )
                    // UnitConversion.g:78:12: VARVALUE
                    {
                    match(input,VARVALUE,FOLLOW_VARVALUE_in_var766); 
                     value = this.varvalue; 

                    }


                    }
                    break;
                case 2 :
                    // UnitConversion.g:79:10: ( VARMW )
                    {
                    // UnitConversion.g:79:10: ( VARMW )
                    // UnitConversion.g:79:12: VARMW
                    {
                    match(input,VARMW,FOLLOW_VARMW_in_var784); 
                     value = this.varmw;    

                    }


                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "var"


    // $ANTLR start "number"
    // UnitConversion.g:83:1: number returns [double value] : NUMBER ;
    public final double number() throws RecognitionException {
        double value = 0.0;

        Token NUMBER17=null;

        try {
            // UnitConversion.g:84:5: ( NUMBER )
            // UnitConversion.g:84:7: NUMBER
            {
            NUMBER17=(Token)match(input,NUMBER,FOLLOW_NUMBER_in_number817); 
             value = Double.parseDouble((NUMBER17!=null?NUMBER17.getText():null)); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }
        finally {
        }
        return value;
    }
    // $ANTLR end "number"

    // Delegated rules


 

    public static final BitSet FOLLOW_sumexpr_in_start43 = new BitSet(new long[]{0x0000000000000000L});
    public static final BitSet FOLLOW_EOF_in_start47 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_mulexpr_in_sumexpr67 = new BitSet(new long[]{0x0000000000000032L});
    public static final BitSet FOLLOW_PLUS_in_sumexpr98 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_mulexpr_in_sumexpr103 = new BitSet(new long[]{0x0000000000000032L});
    public static final BitSet FOLLOW_MINUS_in_sumexpr124 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_mulexpr_in_sumexpr128 = new BitSet(new long[]{0x0000000000000032L});
    public static final BitSet FOLLOW_unaryexpr_in_mulexpr166 = new BitSet(new long[]{0x00000000000000C2L});
    public static final BitSet FOLLOW_MUL_in_mulexpr197 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_unaryexpr_in_mulexpr202 = new BitSet(new long[]{0x00000000000000C2L});
    public static final BitSet FOLLOW_DIV_in_mulexpr223 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_unaryexpr_in_mulexpr228 = new BitSet(new long[]{0x00000000000000C2L});
    public static final BitSet FOLLOW_MINUS_in_unaryexpr279 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_primexpr_in_unaryexpr281 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_primexpr_in_unaryexpr299 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_number_in_primexpr335 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_var_in_primexpr355 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_constant_in_primexpr378 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_funexpr_in_primexpr396 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_parexpr_in_primexpr415 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_E_in_constant455 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_PI_in_constant473 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NA_in_constant490 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_COS_in_funexpr524 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr527 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_SIN_in_funexpr544 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr547 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_SQRT_in_funexpr564 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr566 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_POW_in_funexpr583 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_LPAR_in_funexpr586 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_sumexpr_in_funexpr590 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_COMMA_in_funexpr592 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_sumexpr_in_funexpr596 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_RPAR_in_funexpr598 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_EXP_in_funexpr615 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr618 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_LOG_in_funexpr635 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_LPAR_in_funexpr638 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_sumexpr_in_funexpr642 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_COMMA_in_funexpr644 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_sumexpr_in_funexpr648 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_RPAR_in_funexpr650 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_LOGD_in_funexpr667 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr669 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_LOGN_in_funexpr686 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr688 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_LOGB_in_funexpr705 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_parexpr_in_funexpr707 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_LPAR_in_parexpr736 = new BitSet(new long[]{0x0000000003FCFF20L});
    public static final BitSet FOLLOW_sumexpr_in_parexpr738 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_RPAR_in_parexpr742 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_VARVALUE_in_var766 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_VARMW_in_var784 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NUMBER_in_number817 = new BitSet(new long[]{0x0000000000000002L});

}