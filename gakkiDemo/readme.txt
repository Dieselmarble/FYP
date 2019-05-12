
Abstract:

 OGWE is a computationally optimized fourth-order based ICA/BSS algorithm for the instantaneous case. It uses the SICA algorithm, a quite similar method to the ME-ICA by Comon (1994) with  lower complexity. Besides it does not have the limitation in the number of sources the JADE suffers from. Other contrast may be also minimized with this algorithm: 

AML, AEML, EML, MaSSFOC, ML, MK, SKSE, ...

------------Some more explanations and references are given bellow


OGWE is the optimized jacobi optimization (OJO) applied basically to the SICA algorithm. The real case is considered here.

The SICA algorithm is an algorithm with similar features than the fourth-order based ICA algorithm presented in the well-known paper "Independent component analysis: A new concept?" by P. Comon included Signal Processing, vol. 36, no. 3, pp. 287-314, April 1994. It major advantage is that of a lower computational burden. For futher information on the SICA you may see 

  J.J. Murillo-Fuentes, F.J. González-Serrano. “A sinusoidal contrast function for the separation of statistically independent sources”, Acepted to be published in the IEEE Trans. on Signal Processing. 

or

  J.J. Murillo-Fuentes, F.J. González-Serrano. “Independent component analysis with sinusoidal fourth-order contrasts”, IEEE Int. Conf. on Acoustics, Speech and Signal Processing (ICASSP'2001). IEEE. Salt Lake City (EEUU). Proc. IEEE Int. Conf. on Acoustics, Speech and Signal Processing. pp. 2785-2788. IEEE Press. 2001.

Besides, the n-dimensional case using Jacobi Optimization may be computationally relaxed whenever we have less than 15 sources. This is the Optimized Jacobi Optimization (OJO). This method may be applied to the weighted estimator (WE) by Zarzoso et al. Besides, little improvement of the WE yields the general-WE (GWE) so that a broad family of contrast 

including SICA, AML, AEML, EML, MaSSFOC, ML, MK, SKSE, ... may 

be minimized. The application of OJO to the GWE is explained in detail in

Juan Jose Murillo-Fuentes, Rafael Boloix and Francisco J. González Serrano, “Inititalized Jacobi Optimization in Independent Component Analysis“. In Proc of  ICA’2003. Nara (Japan), Apr 2003. http://viento.us.es/~murillo/ica/ica0.htm

For any further information you may contact me at 
murillo@esi.us.es
http://viento.us.es/~murillo
