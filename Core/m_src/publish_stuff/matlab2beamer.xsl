<?xml version="1.0" encoding="utf-8"?>

<!--
Template for Beamer publish. Creates a simple set of slides. First section  

-->

<!DOCTYPE xsl:stylesheet [ <!ENTITY nbsp "&#160;"> ]>
<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:escape="http://www.mathworks.com/namespace/latex/escape"
  xmlns:mwsh="http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd">
  <xsl:output method="text" indent="no"/>

<xsl:template match="mscript">
    
\documentclass[hyperref, xcolor=dvipsnames, 8pt ]{beamer}
\usepackage{graphbox}    
\usetheme{default}
\usefonttheme[onlymath]{serif}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{items}[default]


\setbeamertemplate{footline}
{
 \includegraphics[align=c, height=0.5cm]{DecyphirLogo.png}%
       \hfill%
       \usebeamercolor[fg]{page number in head/foot}%
       \usebeamerfont{page number in head/foot}%
       \insertframenumber\,/\,\inserttotalframenumber\kern1em%
 }
\usecolortheme[named=BlueViolet]{structure}
\setbeamercovered{transparent=0}

\usepackage[english]{babel}
\usepackage{calc}
\usepackage{amssymb}
\usepackage{amstext}
\usepackage{amsmath}
\usepackage{times}
\usepackage{verbatim}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage{color}
\usepackage{listings}
\usepackage{pgfgantt}
\usepackage[framed]{matlab-prettifier}

\definecolor{lightgray}{gray}{0.5}
\definecolor{bggray}{gray}{0.95}
\lstset{style=Matlab-editor,basicstyle=\scriptsize\mlttfamily, backgroundcolor=\color{bggray}}
\setbeamercolor{code}{fg=black,bg=black!10}
\def \M {\mathcal M}

\begin{document}
    
    <!-- Determine if the there should be an introduction section. -->
    <xsl:variable name="hasIntro" select="count(cell[@style = 'overview'])"/>
    <xsl:if test = "$hasIntro">
\title{<xsl:apply-templates select="cell[1]/steptitle"/>}
\begin{frame}{}
\titlepage

<xsl:apply-templates select="cell[1]/text"/>
</xsl:if>
    
    <xsl:variable name="body-cells" select="cell[not(@style = 'overview')]"/>

    
    <!-- Loop over each cell -->
    <xsl:for-each select="$body-cells">
        <!-- Title of cell -->
        <xsl:if test="steptitle">
          <xsl:variable name="headinglevel">
            <xsl:choose>
              <xsl:when test="steptitle[@style = 'document']">section</xsl:when>
              <xsl:otherwise>begin{frame}[fragile]</xsl:otherwise>
            </xsl:choose>
          </xsl:variable>

\end{frame}
\<xsl:value-of select="$headinglevel"/>{<xsl:apply-templates select="steptitle"/>}

</xsl:if>

        <!-- Contents of each cell -->
        <xsl:apply-templates select="text"/>
        <xsl:apply-templates select="mcode"/>
        <xsl:apply-templates select="mcodeoutput"/>
        <xsl:apply-templates select="img"/>

    </xsl:for-each>


<xsl:if test="copyright">
\footnotesize \color{lightgray} \begin{flushright}
\emph{<xsl:apply-templates select="copyright"/>}
\end{flushright} \color{black} \normalsize
</xsl:if>
\end{frame}
\end{document}

</xsl:template>



<xsl:template name="contents">
  <xsl:param name="body-cells"/>
\subsection{Contents}

\begin{itemize}
\setlength{\itemsep}{-1ex}<xsl:for-each select="$body-cells">
      <xsl:if test="./steptitle">
   \item <xsl:apply-templates select="steptitle"/>
      </xsl:if>
    </xsl:for-each>
\end{itemize}
</xsl:template>




<!-- HTML Tags in text sections -->
<xsl:template match="p">
<xsl:apply-templates/><xsl:text>

</xsl:text>
</xsl:template>

<xsl:template match="ul">\begin{itemize}
<xsl:apply-templates/>\end{itemize}
</xsl:template>
<xsl:template match="ol">\begin{enumerate}
<xsl:apply-templates/>\end{enumerate}
</xsl:template>
<xsl:template match="li">   \item <xsl:apply-templates/><xsl:text>
</xsl:text></xsl:template>

<xsl:template match="pre">
  <xsl:choose>
    <xsl:when test="@class='error'">
\begin{lstlisting}<xsl:value-of select="."/>\end{lstlisting}
    </xsl:when>
    <xsl:otherwise>
\begin{lstlisting}<xsl:value-of select="."/>\end{lstlisting}
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>
<xsl:template match="b">\textbf{<xsl:apply-templates/>}</xsl:template>
<xsl:template match="tt">\texttt{<xsl:apply-templates/>}</xsl:template>
<xsl:template match="i">\textit{<xsl:apply-templates/>}</xsl:template>
<xsl:template match="a">\begin{lstlisting}<xsl:value-of select="."/>\end{lstlisting}</xsl:template>

<xsl:template match="text()">
  <!-- Escape special characters in text -->
  <xsl:call-template name="replace">
    <xsl:with-param name="string" select="."/>
  </xsl:call-template>
</xsl:template>

<xsl:template match="equation">
<xsl:value-of select="."/>
</xsl:template>

<xsl:template match="latex">
    <xsl:value-of select="@text" disable-output-escaping="yes"/>
</xsl:template>

<xsl:template match="html"/>


<!-- Code input and output -->

<xsl:template match="mcode">\begin{lstlisting}
<xsl:value-of select="."/>
\end{lstlisting}
</xsl:template>


<xsl:template match="mcodeoutput">
  <xsl:choose>
    <xsl:when test="substring(.,0,8)='&lt;latex&gt;'">
      <xsl:value-of select="substring(.,8,string-length(.)-16)" disable-output-escaping="yes"/>
    </xsl:when>
    <xsl:otherwise>
      \pause 
      \color{lightgray}
\begin{lstlisting}
<xsl:value-of select="."/>
\end{lstlisting}
\color{black}
</xsl:otherwise>
</xsl:choose>
</xsl:template>

<!-- Figure and model snapshots -->

<xsl:template match="img">
\begin{center}
\visible&lt;2-&gt;{\includegraphics [width=3in]{<xsl:value-of select="@src"/>}}
\end{center}
</xsl:template>

<!-- Colors for syntax-highlighted input code -->

<xsl:template match="mwsh:code">\begin{lstlisting}<xsl:apply-templates/>\end{lstlisting}
</xsl:template>
<xsl:template match="mwsh:keywords">
  <span class="keyword"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:strings">
  <span class="string"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:comments">
  <span class="comment"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:unterminated_strings">
  <span class="untermstring"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:system_commands">
  <span class="syscmd"><xsl:value-of select="."/></span>
</xsl:template>


<!-- Used to escape special characters in the LaTeX output. -->

<escape:replacements>
  <!-- special TeX characters -->
  <replace><from>$</from><to>\$</to></replace>
  <replace><from>&amp;</from><to>\&amp;</to></replace>
  <replace><from>%</from><to>\%</to></replace>
  <replace><from>#</from><to>\#</to></replace>
  <replace><from>_</from><to>\_</to></replace>
  <replace><from>{</from><to>\{</to></replace>
  <replace><from>}</from><to>\}</to></replace>
  <!-- mainly in code -->
  <replace><from>~</from><to>\ensuremath{\tilde{\;}}</to></replace>
  <replace><from>^</from><to>\^{}</to></replace>
  <replace><from>\</from><to>\ensuremath{\backslash}</to></replace>
  <!-- mainly in math -->
  <replace><from>|</from><to>\ensuremath{|}</to></replace>
  <replace><from>&lt;</from><to>\ensuremath{&lt;}</to></replace>
  <replace><from>&gt;</from><to>\ensuremath{&gt;}</to></replace>
</escape:replacements>

<xsl:variable name="replacements" select="document('')/xsl:stylesheet/escape:replacements/replace"/>

<xsl:template name="replace">
  <xsl:param name="string"/>
  <xsl:param name="next" select="1"/>

  <xsl:variable name="count" select="count($replacements)"/>
  <xsl:variable name="first" select="$replacements[$next]"/>
  <xsl:choose>
    <xsl:when test="$next > $count">
      <xsl:value-of select="$string"/>
    </xsl:when>
    <xsl:when test="contains($string, $first/from)">      
      <xsl:call-template name="replace">
        <xsl:with-param name="string"
                        select="substring-before($string, $first/from)"/>
        <xsl:with-param name="next" select="$next+1" />
      </xsl:call-template>
      <xsl:copy-of select="$first/to" />
      <xsl:call-template name="replace">
        <xsl:with-param name="string"
                        select="substring-after($string, $first/from)"/>
        <xsl:with-param name="next" select="$next"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="replace">
        <xsl:with-param name="string" select="$string"/>
        <xsl:with-param name="next" select="$next+1"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
