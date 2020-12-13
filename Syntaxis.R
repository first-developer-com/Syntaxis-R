https://rstudio.com/resources/cheatsheets/
anna.altova@natur.cuni.cz

# BALICKY

intall.packages("package_name") #instalovat balicek
library(package_name) #nacist balicek



# PROMENNE

.x # skryta promenna

# DATOVE OBJEKTY
	# datove typu
	numeric # ciselna hodnota
		integer # a<-5L (je podmnozina numeric)
		double # a<-5.1 (je podmnozina numeric)
	complex # c<-2i
	charakter # retezec
	logical # bool
	na # chybejici hodnota
	NULL # prazdna hodnota
	vector # vektor 
	matrix # matice
	data.frame # data frame
	list # list
	factor # faktor
	casova rada # nebo time series

	#vektor
		#atomicky vektor (ma vzdy hodnoty stejneho typu)
		vektor <- c(1,2,3,4)

	# matice
	m <- matrix(vektor, ncol=x, nrow=x) # matice x stloupce, y radku
	# vektor bude se opkovat pokud matice nebude plna ale pocita se delka vektoru
	# byrow=T plnit matice po radkach (default po stloupcich)
		

class(x) # jaky je datovy typ
is.type # type = charakter/numeric/na (null)/ overeni typu
as.type(x) # prevest k typu as.numeric("5") bude 5


# KONVERTACE DAT


# GRAFY
hist() # histogram
	# parametry hist
		# barva col=""
		# pocet stloupce brakes = int
		# nadpis main = ""

plot(mgp~hp, data=mtcars) # zobraz zavislost stloupce mgp do hp z dataframu mtcars
plot(y=mtcars$mgp, x=mtcars$hp) # zobraz zavislost stloupce mgp do hp z dataframu mtcars
plot(1:20) # nebo tak
	#parametry plot

		# data  x=,y= nebo  x~y, data=DF
		# hlavni nadpis main='' 
		# barva nadpisu col.main = "red"
		# podnadpis sub=''
		# barva podnadpisu col.sub = "yellow"
		# popis osy X xlab = ''
		# popis osy Y ylab = ''
		# rozsah osy X xlim=c(5,15)
		# rozsah osy Y ylim=c(10,15)
		# typy type = 
			# l = carka
			# p = 
			# b = kombinace linie  bodiku
			# c = krouzek
			# h = schody
			# s = 
			# n = 
		# typ bodu pch = inteeger
		# tucnost cary lwd = float
		# typ cary lty = int 
		# barva col="blue" nebo col=1 colors() ukaze zakladni barvy
		# barva cisel na osy col.axis=""
		# barva raamecku fg=""

grid() # pridat daty do grafu ???
lines(2:30) # pridat vrstvu k grafu ???

par(mfrow=c(2,2), mar=c(bot, left, top, rig)) # 2 radky a 1 sploupec tak zobrazovat grafy

	#mfrow matice grafu na obrazky
	# mar margin mezi grafama default (5.1, 4.1, 4.1, 2.1)



barplot() # grafy po stlopceh

data_dp <- c(1, 10, 4, 8, 5)
names(data_dp) <- c("A", "B", "C", "D", "E")
barplot(data_dp, main="GRAF", col = c("blue", "yellow", "green"), legend.text = c("1", "2", "3", "4", "5"))


boxplot(1:100)


# SKRATKY

 CMD + L # Promazat enviroment 



# FUNKCE

help(function) # function
?function # napoveda function
help.search("string") # hledat v napovede string
args(funktion) # zobrazit zakladne argumenty funkce
example(function) #zobrazit vzorec funkce
print(var) # vypsat hodnotu premenny
paste(a,b) # konkatenace s mezerou
paste0(a,b) # konkatenace bez mezery
substr(a, b, c) # orezani radku a=string, b=start, c=to 
letters(c(a:n)) # sequence malych pisem a az n
LETTERS(c(a:n)) # sequence velkych pisem a az n
rm() # vymazat rm(list=ls()) vymazat prostredi
tail() # posledni hodnota promenny
head() # prvni hodnota promenny
system.time(func) # kolil casu pracue funkce
replicate(n,func) # opakovat funk n krat


# IMPORT/EXPORT
getwd() # v jakym adresare jsem
setwd('path/to/dir') # nastavit adresar
dir() # zobraz obsah slozky (argument cesta, kdyz neni tak aktualni)
dir.create("name") # vytvorit adresar
unlink("folder", recurzive=T) # smazat slozku a soubory
file.create("name.ext") # vytvorit soubor
file.remove("name.ext") # smazat soubor


# OPERTORY

	# prirazeni 
		a = "string"  # operator prirazeni nebo prirazeni argumentu funkce
		a <- "string"
	#math	
	+ -	* /
	^ nebo ** # mocnina
	%/% #celociselne deleni
	%% # zbytek deleni (modulo)
	sqrt() # odmocnina
	#logicke
	==, !=, <, >, >=, <=, &, |, !
	c("a", "z") %in% c("a", "b") # element vektoru je v jinym vektoru (T,F)

	

# VEKTOR FUNC

	# operace na vektoru
	* # (1,2,3) * 2 = (2,4,6)
	+ # (1,2,3) + 2 = (3,4,5)
	^ # (1,2,3) ^ 2 = (1,4,9)


a <- 1:10 # [1,2,3,4,5,6,7,8,9,10] generovat vektor 1 az 10
	str(x) # vypsat strukturu objektu
	a[3] # vybrat 3y prvek verktoru
	a[a>5] # vypsat prvky, ktere >5 6,7,8,9,10
	a[a%% == 0] # 2,4,6,8,10 sudi hodnoty
	a[1:3] # hodnoty #1,2,3
	a[c(1,3,5)] # hodnoty #1,3,5
	tail(x,y) # vratit z vektoru x poslednich y hodnot
	length() # kolik prvku ma vektor
	rev(vektor) # otocit vektor


# MATRIX FUNC

dim(m) # dimenze matice (rozmer)
str(m) # do radku
m[x, y] # vybrat x radek, y stloupec
m[1:2, 2:3] # udelat submatice
m[,2] # vybrat druhy stloupec
m[,"second"] # vybrat stloupec second
m[, -c(2)] # matice bez druheho stloupce
attributes(m) # atributy matice 
attributes(m)$dim # atributy dimenze matice 

# kdyz mam vektor 
v <- 1:9
# a udelam 
attributes(v)$dim <- c(3,3)
# tak dostamu mative 3*3

colnames(m) # vypsat jmena stloupce

colnames(m) <- c('first', 'second', 'third')

rownames(m) # vypsat jmena radku

ncol() # pocet stloupce
nrow() # pocet radku
length() # pocet prvku

#maticove operace

A*B # nasobeni element * element
A %*% B # klasicke matematicke nasobeni

rbind(vec_1, vec_2) # spojit matice (prilepit radek)
cbind(vec_1, vec_2) # spojit matice (prilepit stloupec)

rosSums(m) # secist sumy radku
colSums(m) # secist sumy stloupce
rowMeans(m) # prumer radku
colMeans(m) # prumer stloupce
mean(m) # prumer matice
median(m) # median matice

A ><== atd  x # vypise T/F matice pro elementy matice
#   A
# 1,2,3
# 4,5,6
# 7,8,9

A > 5
# F,F,F
# F,F,T
# T,T,T

solve(m) # inverzni matice
det(m) # determinant
t(m) # transponovana matice


# DATA FRAME FUNCTIONS

head(DF, y) # prvni y radku dataframu x
rownames(DF) # zobraz nazvy radku
colnames(DF) # zobraz nazvy stloupce
summary(DF) # zobraz daty o DF
str(DF) # struktura
DF$colname # vypsat jenom stloupec colname z dadaframu DF
DF[3,4] # hodnota 3 radek, 4 stloupec
DF[1,"mpg"] # hodnota 1 radek, stloupec mpg

# vytvorit DF
beatles<-data.frame(
  names=c("John", "Paul", "George", "Ringo"),
  instrument=c("guitar", "bass", "guitar", "drums"),
  birth=c(1940,1942,1943,1940),
  alive=c(F,T,F,T),
  death=c(1980,NA,1943+60,NA)
)

DF$new <- c("") # pridat stloupec z hodnotama

beatles$age_death <- beatles$death - beatles$birth #  vytvorit stloupec age_death a do kazdeho radku vlozit death - birth


nrow(DF) #cislo radku
ncol(DF) # cislo stloupce
rownames(DF) # nazvy radku
colnames(DF) # nazvy stloupce
dim(DF) # dimanze

beatles$age_death[which(beatles$names=="John")] # takovy select filter 
names(DF) # vypsat jmena objektu

merge(X,Y) # spojit 2 DF podle spolecniho stloupce
	# nebo 
	merge(x=X, y=Y, by.x="colname1", by.y="colname2")
	# kdyz stloupce se jmenuji ruzne
	# argument all=TRUE spoji vsechno mozne i kdyz DF1 nema tu polozku v DF2
	# all se da pouzit all.x || all.y
intersect(names(X), names(Y)) # vyhledat spolecne stloupce


# LIST FUNCTIONS

# list muze obsahovat ruzne typy dat(vektoer, tabulka, DF)
# list nema dimenzi

my_list <- list(
  mtcars=mtcars,
  vektor=1:10,
  titanic=Titanic
) #vytvorit list

str(list) # vypsat strukturu listu
length(list) # delka listu

my_list$mtcars # vyber hodnot
my_list$vektor[1] # vyber hodnot

my_list[n] # vrati element pod indexem n jako list

# POZOR my_list[n][n][n][n][n] nikdy nedojde

my_list[[n]] # vrati primo ten element
my_list[["name"]] == my_list$name == my_list[[n]]

names(my_list) # zobrazit/editovat klice listu

# FACTOR FUNCTIONS

faktor <- factor(x=c(1,2,3,4,1,2,4,2), labels=c("prvni", "druhy", "treti", "ctvrty"), ordered = T)
# zalozit vektor
faktor[faktor > "druhy"] # vyselektovat podle podminky

	#netrideny faktor
	barvy <- factor(x=c(2,4,3,2,2,1,1),labels=c("hneda", "zlata", "zluta", "zelena"))

attributes(faktor)
levels(faktor)
as.numeric(faktor) # vypise ten x vektor

# CASOVA RADA FUNKTIONS

# ma parametry
data # daty
start # zacatek v rokach 
end # konec 
frequence # delka periodu 1 = rok, 12 = mesic


# do listu, DF
# vybrat z rady nekrere roky

plot(mesice) # graf
frequency(AirPassengers) # zobraz frequency
start(AirPassengers) # zobraz start
end(AirPassengers) # zobraz end
summary(AirPassengers) # summary
cycle(AirPassengers) # zobraz paradi cyklu
decompose(AirPassengers) # roztridit na vlastnosti
str(decompose(AirPassengers)) 
plot(decompose(AirPassengers)$trend) # grafy vlastnosti
AirPassengers>200





# MATH FUNKTIONS

factorial() #factorial(4) = 4*3*2*1
exp(x) # e^x

cos() # cosinus
sin() # sinus

rnorm() # generovat nahodne cisla normalniho rozdeleni
set.seed() # komplekty nahodnich cisel abych mit stejne nahodne cisla jako nekdo jiny

seq(from=x, to=y, by=z) # generovat sequence on x do y z krokem z
# arg mozny length.out
rep(x,times=y) # opakovat x, y krat args(length.out - opakuj pokud ne vusledek ne bude out, each - kazdu each krat)

mean()
median()
min()
max()
quantile()



# import dat CSV/TXT

carka <- read.table('carka.csv', sep=',', dec='.', header = T, encoding = "UTF-8") # nacita do DF csv soubor sep = oddelavac, dec = desetinny oddelovac, encoding = kodovani

# oddelovac , je "americky format", "europsky format" potrebuje strednik ; 

carka2 <- read.csv('carka.csv') # export csv pro americky formant
strednik2 <- read.csv2('strednik.csv') # export csv pro europsky formant

dalsi moznosti nacitni 
rozdily <- read.table('rozdily.txt', sep='\t', header=T)
rozdily2 <- read.delim('rozdily.txt')


# import dat EXEL

library(readxl) # na to se pouziva balicek
exel <- read_excel('excel.xls'  skip = x) # nahrat tabolku exel, ignorovat prvni x radku

class(exel) # [1] "tbl_df"     "tbl"        "data.frame"

exel<- as.data.frame(exel) # prevest do DF

exel3 <- read_excel('excel2.xlsx', skip=2, sheet='List2', 
                    range='A1:C12', col_names = T, 
                    col_types = c("text", "numeric", "numeric")) # rozsirena varianta

# ukladani grafu

jpeg('rplot.jpeg', width = 500, height = 350, quality = 100) # vytvorit obrazek
plot(prispevek~datum, data = nakaza, type='l') # ulozit graf
dev.off() # ukoncit tvorbu obrazku

png() # analogicke jpeg
 


# export dat z R


write.csv(x, file='soubor.csv', row.names = F, sep='y') # x=data.frame, file = cesta a nazev, row.names = nazvy radku, sep=oddelovac
write.table() # analogicky

save(exel2, file = 'exel2.Rdata') # ulozit
load("/Users/romankrizak/study/Analyza dat a vizualizace/cviceni/2-nd(7.11)/exel2.Rdata") # nacist
dir.create() # vytvorit slozku
save.image(X) # ulozit envirement do souboru
load(x) # nacist envirenemt






# KODOVANI
library(rvest) # knihovna pro kodovani
guess_encoding(x) # jake kodovani ma daty
repair_encoding(x) # NEVIM
Encoding(x) <- "UTF-8" # nastavit kodovani



# rodina APPLY funkce

	# apply() pracuje z pole, matice DF
	# lapply() pracuje z list, vektor. vystupem je vzdy list
	# sapply() vraci vektor, matice, pole. zavisi na situace
	# vapply() verze sapply, ktera umoznuje definovat typ, ktery vraci
	# mapply() cyklicle provadet funkce a jako parameter brat postupne z vektoru
	# replikate() - opakovani 
	# tapply() - provadi funkce pres faktor  kdyz mame treba iris. dataset z tridenim do par skupin
	# by() - jako tapply ale vystup bude typu by
	# aggregate() - podobne jako tapply, vysledkem bude DFnebo casova rada, vstup cokoliv
# to je foreach vektoru

apply(x, margin, fun) # x = vektor, margin = hrana/dimenze, fun = nejaka funkce

dim(mtcars) # [1] 32 11
apply(mtcars,2,min) # na souboru mtcars, pro vsechne stloupce (mar=2), provest funkce min
apply(mtcars,1,mean) # na souboru mtcars, pro vsechne radky (mar=1), ziskat prumer vsech hodnot
apply(mtcars,2,function(x){return (sum(x)*2)}) # apply s anonymni funkce

	# apply na transformace dat
	M <- matrix(1:9,3)
	apply(M,1:2, function(x)x+10)


	# 3-d prostor
	pole <- array(1:27, dim=c(3,3,3))
	apply(pole, c(1,3), sum)

	# argumenty funcke v aaply
	matice = matrix(1:9,3)
	matice[2,2] <- NA
	apply(matice, 2, sum, na.rm=T)

	# na graf
	VADeaths
	par(mfrow=c(2,2))
	apply(VADeaths, 2, barplot)


lapply(list,fun) # list je list nebo vektor, fun je funkce
	# lapply z vektorem
		b<-list(1:10)
		lapply(b[1], mean)

a<-unslist(lapply(list,fun))



# kombinace aaply a lapply
L2 <- list(a=matrix(1:9,3) , b=matrix(10:18,3) , c=matrix(1:9,3))

lapply(L2,FUN=apply,2,sum)


supply(list, fun)
vapply(list, fun, FUN.VAlUE = type) # type jako numeric(1)


Q <- matrix(data=c(rep(1,4), rep(2,4), rep(3,4), rep(4,4)), nrow = 4, ncol = 4) # je to jako priste
Q2 <- mapply(rep, 1:4, 4) # je to jako minule

mapply(function(x,y) x^y, x=c(1,2,3,4), y=c(1,2,3,4)) # 1   4  27 256 1^1,
mapply(function(x) x^2, x=c(1:10)) 1^2, 2^2, 3^2 ..... 


replicate(n=5,expr = sample(1:100,3)) # opakovat 5 krat generovat 3 cisla z 1-100
replicate(n=3, expr=rnorm(10))

tapply(X = iris$Sepal.Length, INDEX = iris$Species, min) # x co pocitam, index = podle ceho tridime, fun = co delame

by(data = mtcars$mpg, INDICES = mtcars$am, mean) # data co pocitam, INDICES = podle ceho tridime, fun = co delame

aggregate(mtcars$mpg, by=list(mtcars$am), mean) # list podle ceho tridit
aggregate(mtcars$mpg, by=list(mtcars$am, mtcars$cyl), mean)

aggregate(AirPassengers, nfrequency = 4, sum) #spocitat sumy passegeru podle kvartalu
aggregate(AirPassengers, nfrequency = 1, sum) #spocitat sumy passegeru podle let
aggregate(AirPassengers, nfrequency = 1, mean) # spocitat prumer passegeru pro  kazdy rok



# cykly
y <- c(6:30)
for (i in y){
  if(i%%7 == 0)
    print(i)
}

x <- TRUE
while(x){
  print('ok')
  x<-FALSE
}


# funkce

my_fun <- function(a,b,c){
  x <<- a+b # dat do globalnich promm hodnotu x
  return(b+c) # vratit
}



#callback
my_fun3 <- function(a){
  a = a+1
  return(a)
}
my_fun4 <- function(b, callback){
  d = mean(b)
  print(d)
  return(callback(d))
}
my_fun4(c(1:15),my_fun3)




# linearni optimalizace

# resi systemy linearnich rovnic
# hlida min/max (maximaliyovat zisk nebo minimalizovat naklady)

1) sestavit linearne rovnice
2) prepsat do matice
3) spocitat matice


Priklad:
  ### Firma KOBLIHY a PNEUMATIKY
  
LIDE
kobliha 1 hodinu
pneumatika 50 hodin
podnikova kapacita 1200 hodin

STROJE
kobliha 2 hodiny
pneumatika 10 hodiny
stoje nesmi bezet dele jak 500 hodin

ZISK
kobliha 10 czk
pneumatika 30 czk
musim vyrobit minimalne 10 ks

CIL: maximalizovat zisk

Rovnice

#  1k  +  50p  <= 1200
#  2k  +  10p  <= 500
#  1k  +  0p   >= 0
#  0k  +  1p   >= 0
#  1k  +  1p   >= 10

To co je vlevo do matice A, co je vpravo, do vektoru B, zisky do vektoru Z

A <- matrix(c(1,50,2,10,1,0,0,1,1,1), nrow = 5, byrow = T)
A
B <- c(1200,500,0,0, 10) 
Z <- c(10,30) # vektor zisku

Pak zadame zmamenka (=, >= ) do vektoru R
R = c('<=', '<=', '>=', '>=', '>=')

library(lpSolve)

vysledek <- lp(
  direction = 'max', # minimalizovat nebo maximalizovat
  objective.in = Z,  # zisk
  const.mat = A,     # naklady
  const.dir = R,     # znamenka > <
  const.rhs = B      # omezeni
  
)
vysledek #vypise max zisk
vysledek$solution #ceho kolik vyrabet



# z semestralky

#Poradi P = Pivo, U = Utopence, K = Klobasa

# Ceny prodejne
Pc = 30 
Uc = 25 
Kc = 40 

# Ceny nakupu
Pnc = 15
Unc = 10
Knc = 12

# Hmotnost
Ph = 0.6
Uh = 0.15
Kh = 0.2

# Vektor Zisku
Z = c(Pc-Pnc, Uc-Unc, Kc-Knc)

# Linearne rovnice
#    1*P +    0*U +    0*K >= 10    Chce mít k dispozici minimálně 10 piv.
#    0*P +    1*U +    1*K >= 10    Chce mít k dispozici minimálně 10 jakýchkoliv jídel.
#    0*P +    1*U +    0*K <= 20    Nesní více než 20 utopenců.
#   Ph*P +   Uh*U +   Kh*K <= 15    Celkem uvezte jen 15 kg potravin!!
#    0*P +    0*U +    1*K <= 30    maximálně 30 klobás na skladě!
#   Pc*P +   Uc*U +   Kc*K <= 2000  O návštěvě víte, že bude disponovat 2.000, - Kč

# Prevedeme do matic pro spracovani lp

A <- matrix(c(1, 0, 0,
              0, 1, 1,
              0, 1, 0,
              Ph,Uh,Kh,
              0, 0, 1,
              Pc,Uc,Kc), ncol = 3, byrow = T)  

B <- c(10, 10, 20, 15, 30, 2000) # hodnoty omezeni

R = c(rep('>=',2),rep('<=',4)) # znamenka

library(lpSolve)

vysledek <- lp(
  direction = 'max', # minimalizovat nebo maximalizovat
  objective.in = Z,  # zisk
  const.mat = A,     # naklady
  const.dir = R,     # znamenka > <
  const.rhs = B      # omezeni
  
)
vysledek #vypise max zisk
vysledek$solution #ceho kolik vyrabet

# Maximalni zisku muze byt 1290 Kc
# Pro to potrebujeme 10 piv, 20 utopence, 30 klobas
# 
# Dam vysledky do rovnici
# 
cbind(t(t(A)*vysledek$solution),R,B)

     
#[1,] "10"  "0"   "0"    ">=" "10"   10 piv      >= 10
#[2,] "0"   "20"  "30"   ">=" "10"   50 jidel    >= 10 
#[3,] "0"   "20"  "0"    "<=" "20"   20 utopence <= 20
#[4,] "6"   "3"   "6"    "<=" "15"   15 kg       <= 15
#[5,] "0"   "0"   "30"   "<=" "30"   30 klobas   <= 30
#[6,] "300" "500" "1200" "<=" "2000" 2000 kc     <= 2000
#
# Odpoved: Pro maximalny zisk 1290 Kc budeme potrebovat prodat 10 piv, 20 utopence, 30 klobas