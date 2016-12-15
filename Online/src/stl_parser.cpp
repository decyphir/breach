// A Bison parser, made by GNU Bison 3.0.4.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// Take the name prefix into account.
#define yylex   CPSGraderlex

// First part of user declarations.
#line 1 "stl_parser.ypp" // lalr1.cc:404

#include "stdafx.h"
#include <iostream>
#include <string>

#include <interval.h>
#include <transducer.h>
#include <signal_expr.h>
#include <robustness.h>
#include <tools.h>

using namespace std;
using namespace CPSGrader;


#line 54 "stl_parser.cpp" // lalr1.cc:404

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

#include "stl_parser.h"

// User implementation prologue.
#line 136 "stl_parser.ypp" // lalr1.cc:412


#include "stl_driver.h"
#include "stl_scanner.h"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex


#line 80 "stl_parser.cpp" // lalr1.cc:412


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (/*CONSTCOND*/ false)
# endif


// Suppress unused-variable warnings by "using" E.
#define YYUSE(E) ((void) (E))

// Enable debugging if requested.
#if CPSGRADERDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << std::endl;                  \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !CPSGRADERDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE(Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void>(0)
# define YY_STACK_PRINT()                static_cast<void>(0)

#endif // !CPSGRADERDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

#line 37 "stl_parser.ypp" // lalr1.cc:479
namespace CPSGrader {
#line 166 "stl_parser.cpp" // lalr1.cc:479

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  Parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr = "";
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              // Fall through.
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }


  /// Build a parser object.
  Parser::Parser (class STLDriver& driver_yyarg)
    :
#if CPSGRADERDEBUG
      yydebug_ (false),
      yycdebug_ (&std::cerr),
#endif
      driver (driver_yyarg)
  {}

  Parser::~Parser ()
  {}


  /*---------------.
  | Symbol types.  |
  `---------------*/

  inline
  Parser::syntax_error::syntax_error (const location_type& l, const std::string& m)
    : std::runtime_error (m)
    , location (l)
  {}

  // basic_symbol.
  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol ()
    : value ()
  {}

  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol (const basic_symbol& other)
    : Base (other)
    , value ()
    , location (other.location)
  {
    value = other.value;
  }


  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const semantic_type& v, const location_type& l)
    : Base (t)
    , value (v)
    , location (l)
  {}


  /// Constructor for valueless symbols.
  template <typename Base>
  inline
  Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const location_type& l)
    : Base (t)
    , value ()
    , location (l)
  {}

  template <typename Base>
  inline
  Parser::basic_symbol<Base>::~basic_symbol ()
  {
    clear ();
  }

  template <typename Base>
  inline
  void
  Parser::basic_symbol<Base>::clear ()
  {
    Base::clear ();
  }

  template <typename Base>
  inline
  bool
  Parser::basic_symbol<Base>::empty () const
  {
    return Base::type_get () == empty_symbol;
  }

  template <typename Base>
  inline
  void
  Parser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move(s);
    value = s.value;
    location = s.location;
  }

  // by_type.
  inline
  Parser::by_type::by_type ()
    : type (empty_symbol)
  {}

  inline
  Parser::by_type::by_type (const by_type& other)
    : type (other.type)
  {}

  inline
  Parser::by_type::by_type (token_type t)
    : type (yytranslate_ (t))
  {}

  inline
  void
  Parser::by_type::clear ()
  {
    type = empty_symbol;
  }

  inline
  void
  Parser::by_type::move (by_type& that)
  {
    type = that.type;
    that.clear ();
  }

  inline
  int
  Parser::by_type::type_get () const
  {
    return type;
  }


  // by_state.
  inline
  Parser::by_state::by_state ()
    : state (empty_state)
  {}

  inline
  Parser::by_state::by_state (const by_state& other)
    : state (other.state)
  {}

  inline
  void
  Parser::by_state::clear ()
  {
    state = empty_state;
  }

  inline
  void
  Parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  inline
  Parser::by_state::by_state (state_type s)
    : state (s)
  {}

  inline
  Parser::symbol_number_type
  Parser::by_state::type_get () const
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[state];
  }

  inline
  Parser::stack_symbol_type::stack_symbol_type ()
  {}


  inline
  Parser::stack_symbol_type::stack_symbol_type (state_type s, symbol_type& that)
    : super_type (s, that.location)
  {
    value = that.value;
    // that is emptied.
    that.type = empty_symbol;
  }

  inline
  Parser::stack_symbol_type&
  Parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    location = that.location;
    return *this;
  }


  template <typename Base>
  inline
  void
  Parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);

    // User destructor.
    switch (yysym.type_get ())
    {
            case 31: // "param_id"

#line 129 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.stringVal); }
#line 419 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 32: // "phi_id"

#line 129 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.stringVal); }
#line 426 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 33: // "new_id"

#line 129 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.stringVal); }
#line 433 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 34: // "signal_id"

#line 129 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.stringVal); }
#line 440 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 35: // "stl_test_id"

#line 129 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.stringVal); }
#line 447 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 36: // "string"

#line 129 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.stringVal); }
#line 454 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 43: // constant_signal

#line 132 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 461 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 44: // signal

#line 131 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 468 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 45: // signal_atom

#line 131 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 475 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 46: // signal_unaryexpr

#line 132 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 482 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 47: // signal_multexpr

#line 132 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 489 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 48: // signal_addexpr

#line 132 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 496 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 49: // signal_expr

#line 132 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 503 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 50: // stl_atom

#line 131 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 510 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 52: // interval

#line 130 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.interval_ptr); }
#line 517 "stl_parser.cpp" // lalr1.cc:617
        break;

      case 53: // stl_formula

#line 131 "stl_parser.ypp" // lalr1.cc:617
        { delete (yysym.value.transducer_ptr); }
#line 524 "stl_parser.cpp" // lalr1.cc:617
        break;


      default:
        break;
    }
  }

#if CPSGRADERDEBUG
  template <typename Base>
  void
  Parser::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " ("
        << yysym.location << ": ";
    YYUSE (yytype);
    yyo << ')';
  }
#endif

  inline
  void
  Parser::yypush_ (const char* m, state_type s, symbol_type& sym)
  {
    stack_symbol_type t (s, sym);
    yypush_ (m, t);
  }

  inline
  void
  Parser::yypush_ (const char* m, stack_symbol_type& s)
  {
    if (m)
      YY_SYMBOL_PRINT (m, s);
    yystack_.push (s);
  }

  inline
  void
  Parser::yypop_ (unsigned int n)
  {
    yystack_.pop (n);
  }

#if CPSGRADERDEBUG
  std::ostream&
  Parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  Parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  Parser::debug_level_type
  Parser::debug_level () const
  {
    return yydebug_;
  }

  void
  Parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // CPSGRADERDEBUG

  inline Parser::state_type
  Parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  inline bool
  Parser::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  inline bool
  Parser::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
  Parser::parse ()
  {
    // State.
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

    // FIXME: This shoud be completely indented.  It is not yet to
    // avoid gratuitous conflicts when merging into the master branch.
    try
      {
    YYCDEBUG << "Starting parse" << std::endl;


    // User initialization code.
    #line 47 "stl_parser.ypp" // lalr1.cc:745
{
    // initialize the initial location object
    yyla.location.begin.filename = yyla.location.end.filename = &driver.streamname;
}

#line 662 "stl_parser.cpp" // lalr1.cc:745

    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, yyla);

    // A new symbol was pushed on the stack.
  yynewstate:
    YYCDEBUG << "Entering state " << yystack_[0].state << std::endl;

    // Accept?
    if (yystack_[0].state == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    // Backup.
  yybackup:

    // Try to take a decision without lookahead.
    yyn = yypact_[yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
        try
          {
            yyla.type = yytranslate_ (yylex (&yyla.value, &yyla.location, driver));
          }
        catch (const syntax_error& yyexc)
          {
            error (yyexc);
            goto yyerrlab1;
          }
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      goto yydefault;

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", yyn, yyla);
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_(yystack_[yylen].state, yyr1_[yyn]);
      /* If YYLEN is nonzero, implement the default value of the
         action: '$$ = $1'.  Otherwise, use the top of the stack.

         Otherwise, the following line sets YYLHS.VALUE to garbage.
         This behavior is undocumented and Bison users should not rely
         upon it.  */
      if (yylen)
        yylhs.value = yystack_[yylen - 1].value;
      else
        yylhs.value = yystack_[0].value;

      // Compute the default @$.
      {
        slice<stack_symbol_type, stack_type> slice (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, slice, yylen);
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
      try
        {
          switch (yyn)
            {
  case 2:
#line 152 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.stringVal) = (yystack_[0].value.stringVal);
        }
#line 774 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 3:
#line 157 "stl_parser.ypp" // lalr1.cc:859
    {
           (yylhs.value.stringVal) = (yystack_[0].value.stringVal);
        }
#line 782 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 4:
#line 163 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.transducer_ptr) = new constant_transducer(*(yystack_[0].value.stringVal));
            (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
            (yylhs.value.transducer_ptr)->param_map = driver.param_map;
            (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;
            delete (yystack_[0].value.stringVal);
        }
#line 794 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 5:
#line 172 "stl_parser.ypp" // lalr1.cc:859
    {
           (yylhs.value.transducer_ptr) = new constant_transducer(*(yystack_[0].value.stringVal));
           (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
           (yylhs.value.transducer_ptr)->param_map = driver.param_map;
           (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;
           delete (yystack_[0].value.stringVal);
        }
#line 806 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 6:
#line 181 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.transducer_ptr) = new signal_transducer(*(yystack_[3].value.stringVal));
            (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
            (yylhs.value.transducer_ptr)->param_map = driver.param_map;
            (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

            int i = driver.signal_map[*(yystack_[3].value.stringVal)];

            if (i==0) {
                cout << "Parsing error: unknown signal " << *(yystack_[3].value.stringVal) << endl;
                delete (yylhs.value.transducer_ptr);
                (yylhs.value.transducer_ptr) = nullptr;
                YYERROR;
            }
            delete (yystack_[3].value.stringVal);
        }
#line 827 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 7:
#line 200 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
        }
#line 835 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 8:
#line 205 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
        }
#line 843 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 9:
#line 210 "stl_parser.ypp" // lalr1.cc:859
    {
	       (yylhs.value.transducer_ptr) = (yystack_[1].value.transducer_ptr);
	    }
#line 851 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 10:
#line 215 "stl_parser.ypp" // lalr1.cc:859
    {
	      (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
	    }
#line 859 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 11:
#line 219 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.transducer_ptr) = new abs_transducer((yystack_[1].value.transducer_ptr));
            (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
            (yylhs.value.transducer_ptr)->param_map = driver.param_map;
            (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

        }
#line 871 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 12:
#line 228 "stl_parser.ypp" // lalr1.cc:859
    {
	      (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
	    }
#line 879 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 13:
#line 232 "stl_parser.ypp" // lalr1.cc:859
    {
	      (yylhs.value.transducer_ptr) = new mult_transducer((yystack_[2].value.transducer_ptr), (yystack_[0].value.transducer_ptr));
          (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
          (yylhs.value.transducer_ptr)->param_map = driver.param_map;
          (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

          }
#line 891 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 14:
#line 241 "stl_parser.ypp" // lalr1.cc:859
    {
	      (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
	    }
#line 899 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 15:
#line 245 "stl_parser.ypp" // lalr1.cc:859
    {
	      (yylhs.value.transducer_ptr) = new plus_transducer((yystack_[2].value.transducer_ptr), (yystack_[0].value.transducer_ptr));
          (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
          (yylhs.value.transducer_ptr)->param_map = driver.param_map;
          (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

          }
#line 911 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 16:
#line 253 "stl_parser.ypp" // lalr1.cc:859
    {
	      (yylhs.value.transducer_ptr) = new minus_transducer((yystack_[2].value.transducer_ptr), (yystack_[0].value.transducer_ptr));
          (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
          (yylhs.value.transducer_ptr)->param_map = driver.param_map;
          (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

          }
#line 923 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 17:
#line 262 "stl_parser.ypp" // lalr1.cc:859
    {
            (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
        }
#line 931 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 18:
#line 268 "stl_parser.ypp" // lalr1.cc:859
    {
              (yylhs.value.transducer_ptr) = new stl_atom((yystack_[2].value.transducer_ptr), *(yystack_[1].value.stringVal), (yystack_[0].value.transducer_ptr));
              (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
              (yylhs.value.transducer_ptr)->param_map = driver.param_map;
              (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;
          }
#line 942 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 19:
#line 276 "stl_parser.ypp" // lalr1.cc:859
    { (yylhs.value.stringVal) = new string("<"); }
#line 948 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 20:
#line 277 "stl_parser.ypp" // lalr1.cc:859
    { (yylhs.value.stringVal) = new string(">"); }
#line 954 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 21:
#line 281 "stl_parser.ypp" // lalr1.cc:859
    {
             (yylhs.value.interval_ptr) = new interval(*(yystack_[3].value.stringVal), *(yystack_[1].value.stringVal));
             delete (yystack_[3].value.stringVal);
             delete (yystack_[1].value.stringVal);
         }
#line 964 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 22:
#line 288 "stl_parser.ypp" // lalr1.cc:859
    {
             (yylhs.value.interval_ptr) = new interval(*(yystack_[2].value.stringVal), *(yystack_[1].value.stringVal));
             delete (yystack_[2].value.stringVal);
             delete (yystack_[1].value.stringVal);
         }
#line 974 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 23:
#line 296 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = (yystack_[0].value.transducer_ptr);
             }
#line 982 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 24:
#line 300 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new not_transducer((yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 994 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 25:
#line 308 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new and_transducer((yystack_[2].value.transducer_ptr), (yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 1006 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 26:
#line 316 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new or_transducer((yystack_[2].value.transducer_ptr), (yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 1018 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 27:
#line 324 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new implies_transducer((yystack_[2].value.transducer_ptr), (yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 1030 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 28:
#line 332 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new ev_transducer((yystack_[1].value.interval_ptr), (yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 1042 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 29:
#line 340 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new alw_transducer((yystack_[1].value.interval_ptr), (yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 1054 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 30:
#line 348 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = new until_transducer((yystack_[3].value.transducer_ptr), (yystack_[1].value.interval_ptr), (yystack_[0].value.transducer_ptr));
                 (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                 (yylhs.value.transducer_ptr)->param_map = driver.param_map;
                 (yylhs.value.transducer_ptr)->signal_map = driver.signal_map;

             }
#line 1066 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 31:
#line 356 "stl_parser.ypp" // lalr1.cc:859
    {
                 (yylhs.value.transducer_ptr) = (yystack_[1].value.transducer_ptr);
             }
#line 1074 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 32:
#line 360 "stl_parser.ypp" // lalr1.cc:859
    {

                 transducer * ref = driver.formula_map[*(yystack_[0].value.stringVal)];

                 if (ref==nullptr) {
                     cout << "Parsing error: unknown identifier " << *(yystack_[0].value.stringVal) << endl;
                     (yylhs.value.transducer_ptr) = nullptr;
                     YYERROR;
                 }
                 else {
                     (yylhs.value.transducer_ptr) = ref->clone();
                     (yylhs.value.transducer_ptr)->trace_data_ptr = &driver.data;
                     delete (yystack_[0].value.stringVal);
                 }
             }
#line 1094 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 33:
#line 379 "stl_parser.ypp" // lalr1.cc:859
    {
                driver.formula_map[*(yystack_[2].value.stringVal)] = (yystack_[0].value.transducer_ptr);
                delete (yystack_[2].value.stringVal);
            }
#line 1103 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 34:
#line 385 "stl_parser.ypp" // lalr1.cc:859
    {
             double x;
             driver.add_trace_test(*(yystack_[2].value.stringVal), *(yystack_[0].value.stringVal), 0., false);
             delete (yystack_[2].value.stringVal);
         }
#line 1113 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 35:
#line 391 "stl_parser.ypp" // lalr1.cc:859
    {
             double x;
             s_to_d(*(yystack_[2].value.stringVal), x);
             driver.add_trace_test(*(yystack_[6].value.stringVal), *(yystack_[4].value.stringVal), x, (yystack_[0].value.boolVal));
             delete (yystack_[6].value.stringVal);
             delete (yystack_[4].value.stringVal);
             delete (yystack_[2].value.stringVal);
         }
#line 1126 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 36:
#line 401 "stl_parser.ypp" // lalr1.cc:859
    {
           (yylhs.value.boolVal) = true;
       }
#line 1134 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 37:
#line 405 "stl_parser.ypp" // lalr1.cc:859
    {
           (yylhs.value.boolVal) = false;
       }
#line 1142 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 38:
#line 412 "stl_parser.ypp" // lalr1.cc:859
    {
                    double val;
                    s_to_d( *(yystack_[0].value.stringVal), val );
                    driver.param_map[*(yystack_[2].value.stringVal)] = val;
                    delete (yystack_[2].value.stringVal);
                    delete (yystack_[0].value.stringVal);
                 }
#line 1154 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 39:
#line 420 "stl_parser.ypp" // lalr1.cc:859
    {
                    double val;
                    s_to_d( *(yystack_[0].value.stringVal), val );
                    driver.param_map[*(yystack_[2].value.stringVal)] = val;
                    delete (yystack_[2].value.stringVal);
                    delete (yystack_[0].value.stringVal);
                 }
#line 1166 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 43:
#line 436 "stl_parser.ypp" // lalr1.cc:859
    {

                    double val;
                    s_to_d( *(yystack_[0].value.stringVal), val );

                    (yylhs.value.param_map_ptr) = new map<string,double>();
                    (*(yylhs.value.param_map_ptr))[*(yystack_[2].value.stringVal)] = val;
                    delete (yystack_[2].value.stringVal);
                    delete (yystack_[0].value.stringVal);
                 }
#line 1181 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 44:
#line 447 "stl_parser.ypp" // lalr1.cc:859
    {
                    double val;
                    s_to_d( *(yystack_[0].value.stringVal), val );
                    (yylhs.value.param_map_ptr) = new map<string,double>();
                    (*(yylhs.value.param_map_ptr))[*(yystack_[2].value.stringVal)] = val;
                    delete (yystack_[2].value.stringVal);
                    delete (yystack_[0].value.stringVal);
                 }
#line 1194 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 46:
#line 459 "stl_parser.ypp" // lalr1.cc:859
    {
    (yylhs.value.param_map_ptr) = new map<string,double>(*(yystack_[2].value.param_map_ptr));
    auto elem = (yystack_[0].value.param_map_ptr)->begin();
    (*(yylhs.value.param_map_ptr))[elem->first] = elem->second;
}
#line 1204 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 47:
#line 466 "stl_parser.ypp" // lalr1.cc:859
    {
                           (yylhs.value.param_map_ptr) = new map<string,double>();
                       }
#line 1212 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 48:
#line 470 "stl_parser.ypp" // lalr1.cc:859
    {
                           (yylhs.value.param_map_ptr) = (yystack_[1].value.param_map_ptr);
                       }
#line 1220 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 49:
#line 477 "stl_parser.ypp" // lalr1.cc:859
    {
             short idx =  driver.signal_map.size()+1;
             driver.signal_map[*(yystack_[0].value.stringVal)] = idx;
             delete (yystack_[0].value.stringVal);
          }
#line 1230 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 50:
#line 483 "stl_parser.ypp" // lalr1.cc:859
    {
              delete (yystack_[0].value.stringVal);
          }
#line 1238 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 54:
#line 493 "stl_parser.ypp" // lalr1.cc:859
    {
            driver.add_stl_test(*(yystack_[10].value.stringVal), *(yystack_[9].value.param_map_ptr), (yystack_[7].value.transducer_ptr), *(yystack_[5].value.stringVal), *(yystack_[3].value.stringVal), (yystack_[1].value.boolVal) );
            delete (yystack_[10].value.stringVal);
            delete (yystack_[9].value.param_map_ptr);
            delete (yystack_[5].value.stringVal);
            delete (yystack_[3].value.stringVal);
        }
#line 1250 "stl_parser.cpp" // lalr1.cc:859
    break;

  case 55:
#line 501 "stl_parser.ypp" // lalr1.cc:859
    {
            if ((yystack_[0].value.param_map_ptr)->empty())
                driver.add_stl_test(*(yystack_[1].value.stringVal));
            else {
                driver.add_stl_test(*(yystack_[1].value.stringVal), *(yystack_[0].value.param_map_ptr));
            }
            delete (yystack_[1].value.stringVal);
            delete (yystack_[0].value.param_map_ptr);
        }
#line 1264 "stl_parser.cpp" // lalr1.cc:859
    break;


#line 1268 "stl_parser.cpp" // lalr1.cc:859
            default:
              break;
            }
        }
      catch (const syntax_error& yyexc)
        {
          error (yyexc);
          YYERROR;
        }
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, yylhs);
    }
    goto yynewstate;

  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yyla.location, yysyntax_error_ (yystack_[0].state, yyla));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;
    yyerror_range[1].location = yystack_[yylen - 1].location;
    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[yystack_[0].state];
          if (!yy_pact_value_is_default_ (yyn))
            {
              yyn += yyterror_;
              if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
                {
                  yyn = yytable_[yyn];
                  if (0 < yyn)
                    break;
                }
            }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yyerror_range[1].location = yystack_[0].location;
          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = yyn;
      yypush_ ("Shifting", error_token);
    }
    goto yynewstate;

    // Accept.
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    // Abort.
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack"
                 << std::endl;
        // Do not try to display the values of the reclaimed symbols,
        // as their printer might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
  }

  void
  Parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what());
  }

  // Generate an error message.
  std::string
  Parser::yysyntax_error_ (state_type yystate, const symbol_type& yyla) const
  {
    // Number of reported tokens (one for the "unexpected", one per
    // "expected").
    size_t yycount = 0;
    // Its maximum.
    enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
    // Arguments of yyformat.
    char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];

    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state
         merging (from LALR or IELR) and default reductions corrupt the
         expected token list.  However, the list is correct for
         canonical LR with one exception: it will still contain any
         token that will not be accepted due to an error action in a
         later state.
    */
    if (!yyla.empty ())
      {
        int yytoken = yyla.type_get ();
        yyarg[yycount++] = yytname_[yytoken];
        int yyn = yypact_[yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            // Stay within bounds of both yycheck and yytname.
            int yychecklim = yylast_ - yyn + 1;
            int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
            for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
              if (yycheck_[yyx + yyn] == yyx && yyx != yyterror_
                  && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
                {
                  if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                    {
                      yycount = 1;
                      break;
                    }
                  else
                    yyarg[yycount++] = yytname_[yyx];
                }
          }
      }

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
        YYCASE_(0, YY_("syntax error"));
        YYCASE_(1, YY_("syntax error, unexpected %s"));
        YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    size_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += yytnamerr_ (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char Parser::yypact_ninf_ = -84;

  const signed char Parser::yytable_ninf_ = -1;

  const signed char
  Parser::yypact_[] =
  {
     -16,    -4,    38,   -20,    41,   -84,    28,   -84,   -84,   -84,
       2,    26,    48,   -84,    75,   -84,   -84,   -84,    83,    54,
      25,    12,   -84,   -84,   -84,   -84,    71,    72,    -4,    38,
      67,    25,    25,   101,   101,   100,   -84,   -84,   -84,   103,
     -84,   -84,   -84,   -84,   -84,    74,    27,   -84,    78,   102,
     102,   -84,   -19,   -84,   -84,   -84,   -84,    91,    45,    64,
      10,    47,    25,    25,    35,    93,    19,    19,    19,   -84,
     -84,    35,    25,    25,    25,   101,    30,    73,   -84,   -84,
     -84,    80,   -84,   -84,   -84,   -84,   -11,   -84,   -84,    35,
     105,   109,   -84,   -84,   -84,   -84,     0,     0,     0,    25,
      76,    77,   -84,     1,    25,    99,    47,   113,   112,   -84,
     -84,   -84,    89,    90,   -84,    30,    70,    56,   117,   -84,
     -84,   -84,   -84,    86,   -84,   -84,   -84,   -84,   106,    88,
     108,    56,    87,   -84
  };

  const unsigned char
  Parser::yydefact_[] =
  {
       0,     0,     0,     0,     0,    59,     0,    61,    62,    60,
       0,     0,     0,    40,    42,    49,    50,    51,    53,     0,
       0,     0,     1,    63,    65,    64,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     4,     5,    32,     0,
       8,     7,    10,    12,    14,    17,     0,    23,    33,    47,
      47,    57,     0,    38,    39,    41,    52,    34,     0,     0,
      24,     0,     0,     0,     0,     0,     0,     0,     0,    19,
      20,     0,     0,     0,     0,     0,     0,     0,    55,    58,
      56,     0,     9,    31,     2,     3,     0,    29,    28,     0,
       0,     0,    15,    16,    13,    18,    25,    26,    27,     0,
       0,     0,    45,     0,     0,     0,     0,     0,     0,    11,
       6,    30,     0,     0,    48,     0,     0,     0,     0,    22,
      43,    44,    46,     0,    36,    37,    35,    21,     0,     0,
       0,     0,     0,    54
  };

  const signed char
  Parser::yypgoto_[] =
  {
     -84,   -83,   -84,   -84,    31,   -84,   -84,   -84,   -27,   -84,
     -84,   -29,   -31,   116,   -84,    -3,   104,   -84,   119,    15,
     -84,    81,   107,   -84,   -84,    82,   -84,   123,   -84
  };

  const signed char
  Parser::yydefgoto_[] =
  {
      -1,    86,    40,    41,    42,    43,    44,    45,    46,    47,
      71,    62,    48,     5,     6,   126,    13,    14,     7,   102,
     103,    77,    17,    18,     8,    51,    52,     9,    10
  };

  const unsigned char
  Parser::yytable_[] =
  {
      59,    60,    22,   107,    58,    63,   106,   114,     1,     2,
       3,    73,    74,    19,    49,    75,    50,     4,   115,    84,
      85,    79,    74,   118,    89,    75,     1,    11,     3,    12,
      31,    87,    88,    32,    69,     4,    70,    90,    33,    34,
      89,    96,    97,    98,    95,    49,    99,    50,    35,    36,
      37,    82,    69,    39,    70,    36,    37,    38,    35,    39,
      20,   100,   108,   101,    26,    36,    37,    21,   111,    39,
      83,    15,    16,   116,    72,    73,    74,    84,    85,    75,
      72,    73,    74,   124,   125,    75,    27,   123,    72,    73,
      74,    30,    28,    75,    66,    67,    68,    92,    93,    94,
      29,    53,    54,    57,    61,    64,    65,    76,    81,    91,
     105,   109,   104,   110,   112,   113,   117,   119,    82,   120,
     121,   127,   128,   129,   130,   131,    23,   133,   132,    24,
     122,    78,    55,    25,    80,     0,    56
  };

  const short int
  Parser::yycheck_[] =
  {
      31,    32,     0,    86,    31,    34,    17,     6,    24,    25,
      26,    11,    12,    33,    33,    15,    35,    33,    17,    30,
      31,    40,    12,   106,     5,    15,    24,    31,    26,    33,
       5,    62,    63,     8,     7,    33,     9,    64,    13,    14,
       5,    72,    73,    74,    71,    33,    75,    35,    23,    30,
      31,     6,     7,    34,     9,    30,    31,    32,    23,    34,
      19,    31,    89,    33,    38,    30,    31,    39,    99,    34,
       6,    33,    34,   104,    10,    11,    12,    30,    31,    15,
      10,    11,    12,    27,    28,    15,    38,    17,    10,    11,
      12,    37,    17,    15,    20,    21,    22,    66,    67,    68,
      17,    30,    30,    36,     3,     5,     3,     5,    17,    16,
      30,     6,    39,     4,    38,    38,    17,     4,     6,    30,
      30,     4,    36,    17,    36,    17,    10,    40,   131,    10,
     115,    50,    28,    10,    52,    -1,    29
  };

  const unsigned char
  Parser::yystos_[] =
  {
       0,    24,    25,    26,    33,    54,    55,    59,    65,    68,
      69,    31,    33,    57,    58,    33,    34,    63,    64,    33,
      19,    39,     0,    54,    59,    68,    38,    38,    17,    17,
      37,     5,     8,    13,    14,    23,    30,    31,    32,    34,
      43,    44,    45,    46,    47,    48,    49,    50,    53,    33,
      35,    66,    67,    30,    30,    57,    63,    36,    49,    53,
      53,     3,    52,    52,     5,     3,    20,    21,    22,     7,
       9,    51,    10,    11,    12,    15,     5,    62,    62,    40,
      66,    17,     6,     6,    30,    31,    42,    53,    53,     5,
      49,    16,    45,    45,    45,    49,    53,    53,    53,    52,
      31,    33,    60,    61,    39,    30,    17,    42,    49,     6,
       4,    53,    38,    38,     6,    17,    53,    17,    42,     4,
      30,    30,    60,    17,    27,    28,    56,     4,    36,    17,
      36,    17,    56,    40
  };

  const unsigned char
  Parser::yyr1_[] =
  {
       0,    41,    42,    42,    43,    43,    44,    45,    45,    45,
      46,    46,    47,    47,    48,    48,    48,    49,    50,    51,
      51,    52,    52,    53,    53,    53,    53,    53,    53,    53,
      53,    53,    53,    54,    55,    55,    56,    56,    57,    57,
      58,    58,    59,    60,    60,    61,    61,    62,    62,    63,
      63,    64,    64,    65,    66,    66,    67,    67,    68,    69,
      69,    69,    69,    69,    69,    69
  };

  const unsigned char
  Parser::yyr2_[] =
  {
       0,     2,     1,     1,     1,     1,     4,     1,     1,     3,
       1,     4,     1,     3,     1,     3,     3,     1,     3,     1,
       1,     5,     4,     1,     2,     3,     3,     3,     3,     3,
       4,     3,     1,     3,     4,     8,     1,     1,     3,     3,
       1,     3,     2,     3,     3,     1,     3,     0,     3,     1,
       1,     1,     3,     2,    11,     2,     2,     1,     4,     1,
       1,     1,     1,     2,     2,     2
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const Parser::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "\"[\"", "\"]\"", "\"(\"",
  "\")\"", "\"<\"", "\"not\"", "\">\"", "\"and\"", "\"or\"", "\"=>\"",
  "\"alw\"", "\"ev\"", "\"until\"", "\"time\"", "\",\"", "\"end of line\"",
  "\":=\"", "\"+\"", "\"-\"", "\"*\"", "\"abs\"", "\"param_decl\"",
  "\"signal_decl\"", "\"test\"", "\"true\"", "\"false\"", "\"double\"",
  "\"constant\"", "\"param_id\"", "\"phi_id\"", "\"new_id\"",
  "\"signal_id\"", "\"stl_test_id\"", "\"string\"", "':'", "'='", "'{'",
  "'}'", "$accept", "constant", "constant_signal", "signal", "signal_atom",
  "signal_unaryexpr", "signal_multexpr", "signal_addexpr", "signal_expr",
  "stl_atom", "op", "interval", "stl_formula", "assignement", "trace_env",
  "boolean", "param_assignement", "param_assignement_list",
  "param_assignements", "local_param_assignement",
  "local_param_assignement_list", "local_param_assignements", "signal_new",
  "signal_new_list", "signal_decl", "stl_test", "stl_tests", "trace_test",
  "start", YY_NULLPTR
  };

#if CPSGRADERDEBUG
  const unsigned short int
  Parser::yyrline_[] =
  {
       0,   151,   151,   156,   162,   171,   180,   199,   204,   209,
     214,   218,   227,   231,   240,   244,   252,   261,   267,   276,
     277,   280,   287,   295,   299,   307,   315,   323,   331,   339,
     347,   355,   359,   378,   384,   390,   400,   404,   411,   419,
     428,   429,   431,   435,   446,   457,   458,   466,   469,   476,
     482,   487,   488,   490,   492,   500,   511,   512,   514,   516,
     517,   518,   519,   520,   521,   522
  };

  // Print the state stack on the debug stream.
  void
  Parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << i->state;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  Parser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):" << std::endl;
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // CPSGRADERDEBUG

  // Symbol number corresponding to token number t.
  inline
  Parser::token_number_type
  Parser::yytranslate_ (int t)
  {
    static
    const token_number_type
    translate_table[] =
    {
     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    37,     2,
       2,    38,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    39,     2,    40,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36
    };
    const unsigned int user_token_number_max_ = 291;
    const token_number_type undef_token_ = 2;

    if (static_cast<int>(t) <= yyeof_)
      return yyeof_;
    else if (static_cast<unsigned int> (t) <= user_token_number_max_)
      return translate_table[t];
    else
      return undef_token_;
  }

#line 37 "stl_parser.ypp" // lalr1.cc:1167
} // CPSGrader
#line 1781 "stl_parser.cpp" // lalr1.cc:1167
#line 523 "stl_parser.ypp" // lalr1.cc:1168


void CPSGrader::Parser::error(const Parser::location_type& l,
			    const std::string& m)
{
    driver.error(l, m);
}
