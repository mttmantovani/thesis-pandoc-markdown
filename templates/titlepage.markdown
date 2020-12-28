\pagenumbering{roman}
\setcounter{page}{1}
\begin{titlepage}
    \begin{center}
        \vspace*{6mm}
        \fontsize{25}{30}\selectfont
        \textbf{$title-alt$}\\[2.0cm]
        \fontsize{11}{13.2}\selectfont
        {\Large \textbf{Dissertation zur Erlangung des} \\
         \textbf{akademischen Grades eines Doktors} \\
         \textbf{der Naturwissenschaften} \\
         \textbf{(Dr.rer.nat.)}} \\[1.5cm]
        {\large vorgelegt von}  \\[11pt]
        {\large $surname$, $name$} \\[1.5cm]
        {\large an der} \\[11pt]
        \includegraphics[height=2cm]{./figures/unilogo.pdf} \\ [2cm]
        {\large Mathematisch-Naturwissenschaftliche Sektion} \\
        {\large Fachbereich Physik}
    \end{center}
    \vfill
    {\large $date$}
\end{titlepage}

$if(print)$
\blankpage
$if(published)$
\thispagestyle{empty}
\vspace*{\fill}
\begin{defenseinfo}
Tag der mündlichen Prüfung: $defense-date$ \\
1. Referent: $first-referee$ \\
2. Referent: $second-referee$ \\
\end{defenseinfo}
$else$
\blankpage
$endif$
\blankpage
$else$
$if(published)$
\clearpage
\thispagestyle{empty}
\vspace*{\fill}
\begin{defenseinfo}
Tag der mündlichen Prüfung: $defense-date$ \\
1. Referent: $first-referee$ \\
2. Referent: $second-referee$ \\
\end{defenseinfo}
\cleardoublepage
$else$
\cleardoublepage
$endif$
$endif$

