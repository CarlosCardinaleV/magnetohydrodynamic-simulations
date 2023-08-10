# Magnetohydrodynamic (MHD) Simulations

# Contexto del codigo
El codigo es una simulación del paper llamado "Implicit predictor–corrector central finite difference scheme for the equations of magnetohydrodynamic simulations" hecho por los siguientes investigadores:
1. T.C. Tsai
2. H.-S. Yu
3. M.-S. Hsieh
4. S.H. Lai
5. Y.-H. Yang

En el abstract del paper explican que el artículo propone un esquema de diferencia finita central implícita de alto orden predictor-corrector (iPCCFD) y demuestra su alta eficiencia en el cómputo paralelo. De especial interés son los estudios numéricos a gran escala, como las simulaciones magnetohidrodinámicas (MHD) en la magnetosfera planetaria. Se desarrolla un esquema iPCCFD basado en el método de diferencia finita central de quinto orden y el método implícito predictor-corrector de cuarto orden en combinación con la técnica de eliminación de errores de redondeo (ERE). Se examinan varios estudios numéricos, como el problema unidimensional del tubo de choque Brio-Wu, el sistema de vórtices Orszag-Tang bidimensional, la inestabilidad de tipo K-H de vórtices, la inestabilidad de tipo K-H de tipo kink, la advección de lazo de campo y la onda de choque. Todos los resultados de las simulaciones son consistentes con numerosas literaturas. iPCCFD puede minimizar las inestabilidades y ruidos numéricos junto con los términos adicionales de difusión. Todos nuestros estudios presentan errores numéricos relativamente pequeños sin emplear ninguna reconstrucción libre de divergencia. En particular, obtenemos resultados bastante estables en el problema bidimensional del tubo de choque Brio-Wu que conserva bien ∇ · B = 0 en toda la simulación. La técnica ERE elimina la acumulación de errores de redondeo en el sistema uniforme o no perturbado. También hemos demostrado que iPCCFD se caracteriza por su alto orden de precisión y baja disipación numérica en las pruebas de onda de Alfvén polarizada circularmente. El esquema propuesto iPCCFD es un esquema numérico paralelo-eficiente y de alta precisión para resolver las ecuaciones MHD en sistemas de conservación hiperbólica.
