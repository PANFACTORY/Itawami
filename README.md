# iTawami
iTawami shows deformation of a rectangle plate of which all edges are supported or fixed with central concentrated force.
When calculating a plate of which all edges are supported, iTawami expresses the deflection w as:

<img src="https://latex.codecogs.com/gif.latex?w\left&space;(&space;x,y&space;\right&space;)=&space;\frac{1}{\pi&space;^{4}D}\sum_{m=1}^{\infty&space;}\sum_{n=1}^{\infty&space;}\frac{a_{mn}}{\left&space;(&space;\left&space;(&space;\frac{m}{a}&space;\right&space;)^{2}&plus;\left&space;(&space;\frac{n}{b}&space;\right&space;)^{2}&space;\right&space;)^{2}}sin\frac{m\pi&space;x}{a}sin\frac{n\pi&space;y}{b}"/>

and 

<img src="https://latex.codecogs.com/gif.latex?a_{mn}=&space;\frac{4}{ab}\int_{0}^{a}\int_{0}^{b}q\left&space;(&space;x,y&space;\right&space;)sin\frac{m\pi&space;x}{a}sin\frac{n\pi&space;y}{b}dxdy"/>

It is known as Navier's solution. When calculating a plate of which all edges are fixed, iTawami assumes the deflection w as:

<img src="https://latex.codecogs.com/gif.latex?w\left&space;(&space;x,y&space;\right&space;)=&space;\sum_{m=1}^{\infty&space;}\sum_{n=1}^{\infty&space;}c_{mn}\left&space;(&space;1-cos\frac{2\pi&space;mx}{a}&space;\right&space;)\left&space;(&space;1-cos\frac{2\pi&space;ny}{b}&space;\right&space;)"/>

and detects c based on Ritz's method.


## Dependency
iTawami runs on Python3 and mainly depends on these python packages:
* math
* numpy
* tkinter
* matplotlib
* scipy
* numba


## Usage
iTawami is operated with pressing these keys:
* Å™Å´Å®Å© : rotate graph
* s : show a plate of which all edges are supported 
* f : show a plate of which all edges are fixed
* a : change width a of plate after typing digit keys
* b : change length b of plate after typing digit keys


## Author
Tanabe Yuta


## References
* S.Timoshenko and S.Woinowsky-Krieger: THEORY OF PLATES AND SHELLS, McGRAW-HILL BOOK COMPANY(1959)
* Hangai Yasuhiko: ïΩî¬ÇÃäÓëbóùò_, è≤çëé–(1995)



