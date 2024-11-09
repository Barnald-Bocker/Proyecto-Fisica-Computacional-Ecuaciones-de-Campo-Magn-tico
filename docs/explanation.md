#Metodos numericos usados en el proyecto 

##Metodo de relajacion de jacobi 

Este es un metodo que utiliza diferencias finitas para definir las derivadas del siguiente problema:

\begin{equation}
   \frac{\partial^{2}\phi}{\partial{x^{2}}} + \frac{\partial^{2}\phi}{\partial{y^{2}}} = 0
\end{equation}

en donde las derivadas vistas se pueden aproximar usando diferencias centrales para asi poderlas expresar de la siguiente manera

\begin{equation}
 \frac{\partial^{2}\phi}{\partial{x^{2}}} \approx \frac{\phi (x +a, y) - 2\phi (x,y) + \phi (x -a,y)}{a^{2}}  
\end{equation}

\begin{equation}
 \frac{\partial^{2}\phi}{\partial{y^{2}}} \approx \frac{\phi (x,y + a) - 2\phi (x,y) + \phi (x,y - a)}{a^{2}} 
\end{equation} 

De esta manera, la ecuacion diferencial dada se puede aproximar de la siguiente manera:
\begin{align}
 &\frac{\phi (x +a, y) - 2\phi (x,y) + \phi (x -a,y)}{a^{2}} + \frac{\phi (x,y + a) - 2\phi (x,y) + \phi (x,y - a)}{a^{2}} \approx  0 \\
 &\frac{\phi (x +a, y) + \phi (x -a,y)}{a^{2}} + \frac{\phi (x,y + a) + \phi (x,y - a)}{a^{2}} -4\phi (x,y) \approx 0 \\
 
\end{align} 
