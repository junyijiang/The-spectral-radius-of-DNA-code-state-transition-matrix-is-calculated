# The spectral radius of DNA code state transition matrix is calculated
This project is a record of the search processes for DNA secondary structure avoidance codes.

Notation:

sequence: DNA string over D={A,C,G,T}. 

                For example x = ACATCTGA.

m_subsequence: sub-sequence x length of m. 

                For example m = 3. a sub-sequence y = ACA.

Rc(x): reverse complementary sequence of x with Watson-Crick(A-T,T-A,C-G,G-C). 
                
                For example Rc(x) = TCAGATGT.

xRcx: a pair of sub-sequences x and y, with y = Rc(x).

replace sequence: {A,T,C,G} ——> {T,A,G,C}. 
                
                For example replace of x: TGTAGACT. (the replace rule not unique)

1.avoid_bases.

Count all m-xRcx. For example, when m = 2, m-xRcx: 

                all x :        AA AC AG CA CC GA
                all Rc(x) :    TT GT CT TG GG TC
                
Select any one of xRcx according the rules. For example: we get a set of m-sub-sequences: {AA, GT, AG, TG, CC, GA}.

                all x :        AA    AG    CC GA
                all Rc(x) :       GT    TG             
            
Construct a state transition matrix based on all selected sequences.

                  AA GT AG TG CC GA
              AA  1  0  1  0  0  0
              GT  0  0  0  1  0  0
              AG  0  1  0  0  0  1
              TG  0  1  0  0  0  1
              CC  0  0  0  0  1  0
              GA  1  0  1  0  0  0

Count the matrix maximum eigenvalue.

              max_eigenvalue： 1.801937735804836072

2.delete_replacement.

We have seven replace rules: {TACG, TAGC, ATGC, CGAT, CGTA, GCAT, GCTA}

For every set of xRcx combinations. cout seven replace set, and delete some one set equal to xRcx combinations.
