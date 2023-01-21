####### WORKSHOP DIA 1: Exercício básico de R ########


#Basics of R
#FUNCTIONS
#functions are things we tell R to do for us, such as arithmetic or printing
#  words.

#Run the following lines
2+2
print("Hello world")

#ARGUMENTS 
# Things we need to supply to R functions for them to run properly
# print("Hello world") <- what is inside the "" is the argument for the function "print"

#OBJECTS

# We want to use R to take the average of 1 2 3 4 and 5

mean(1,2,3,4,5) #<- we want to take the average of this series, but R can't read
#past the first number

# We need to put lists of numbers in vector format c() for R to recognize them
mean(c(1,2,3,4,5))

#But what if our data was 1000s of entries long?
#Instead of typing it out, we will assign the list to an object called "mynumbers"
mynumbers <- c(1,2,3,4,5) # <- is the assigner
mynumbers #enter this, what does R output?
mean(mynumbers) #Now we can take the same mean using this new object instead of
#typing out the whole vector

#Ways to store numeric values

#1. Vectors (linear list, 1 dimensional)
mynumbers<- c(1,2,3,4,5)

#2. Data Frames ( x,y list, 2 dimensional)
first_column <- c("value_1", "value_2", "value_3")
second_column <- c("value_3", "value_4", "value_5")

mydataframe <- data.frame(first_column, second_column) #data.frame is the function
mydataframe

#Object classes

# In R, we can work with various types of data (numbers, names, variables, etc)
# How does R know the difference?

#1. Numerics
1,2,3,4,5
class(1) #<-class() function tells us what category the number 1 belongs to
class("1") #<- by enclosing numbers or names in quotes, R recognizes them as characters
class(mydataframe) #R recognizes this object as a data frame!

# EXAMPLES
mean(c(1,2,3,4,5))
#[1] 3

#2. CHARACTERS
mean(c("1","2","3","4","5"))
#[1] NA
#  Warning message:
#  In mean.default(c("1", "2", "3", "4", "5")) :
#  argument is not numeric or logical: returning NA    
# R cannot perform the mean on things that are not numbers!

#3 LOGICALS
#TRUE, FALSE
#These are basically matching expressions where R checks if the result 
#is true or false. Try the expression below:

#Try the following lines of code. What do you get?
2==3


2==2

#########Subsetting dataframes##########
#For this we will use the pre-loaded dataframe, mtcars

data(mtcars)

#mtcars has a lot of info, we don't need it all. Let's tell R to 
#subset the variables "mpg" and "hp".

#Here is how to index (subset) rows and columns of a dataframe
#df[rows,columns]

#Let's do mpg first
mpg<-mtcars[,"mpg"]
mpg<-mtcars[,1]
#we have the list of numbers, but we are missing the information about the cars
#they came from. How can we grab those?

#Here are the names from our dataframe
rownames(mtcars)

#The function below will assign them to the mpg vector
#names() assigns names to VECTORS
#rownames() pulls names of rows in DATAFRAME
names(mpg)<-rownames(mtcars)
mpg

#######We want info on both hp, wt, and mpg, but we don't want them as two separate
#vectors. Let's create a dataframe using the subsetting we learned above

#run these and see what happens. Why can R only correctly understand the second one?
correlation.df<-mtcars[,"mpg", "hp", "wt"]
correlation.df<- mtcars[,c("mpg", "hp", "wt")]


############### Now we can plot a regression to test the hypothesis that
###### car weight or hp predicts mpg!

correlation.df<- mtcars[,c("mpg", "hp", "wt")]

#Structure of a regression;
#y = ax + b
#Following is the description of the parameters used −
#y is the response variable.

#x is the predictor variable.

#a and b are constants which are called the coefficients.

x <- correlation.df$wt
y <- correlation.df$mpg
z<- correlation.df$hp

# Plot with main and axis titles
# Change point shape (pch = 19) and remove frame.
plot(x, y, main = "Relationship between mpg and wt",
     xlab = "Vehicle Weight", ylab = "MPG",
     pch = 19, frame = FALSE)
# Add regression line
plot(x, y, main = "Relationship between mpg and wt",
     xlab = "Vehicle Weight", ylab = "MPG",
     pch = 19, frame = FALSE)
abline(lm(y ~ x, data = correlation.df), col = "blue")

#abline adds a regression line to the plot
# lm = function for linear model
# lm is specified as lm(y ~ x, data = dataframe)

#Maybe we would like to use our plot outside of RStudio. We can write it to a PDF!
#let's save our beautiful plot using the PDF function
pdf(file = "WeightRegression.pdf")
plot(x, y, main = "Relationship between mpg and wt",
     xlab = "Vehicle Weight", ylab = "MPG",
     pch = 19, frame = FALSE)
abline(lm(y ~ x, data = correlation.df), col = "blue")
dev.off()

#How can we retrieve the statistical information about our plot?
lm(y ~ x, data = correlation.df)

# Call:
#lm(formula = y ~ x, data = correlation.df)

#Coefficients:
#  (Intercept)            x  
#   37.285             -5.344 


#Now let's plot the relationship of hp and mpg
# Plot with main and axis titles
# Change point shape (pch = 19) and remove frame.
plot(z, y, main = "Relationship between mpg and hp",
     xlab = "Vehicle Weight", ylab = "HP",
     pch = 19, frame = FALSE)
# Add regression line
plot(z, y, main = "Relationship between mpg and hp",
     xlab = "Vehicle Weight", ylab = "HP",
     pch = 19, frame = FALSE)
abline(lm(y ~ z, data = correlation.df), col = "red")

lm(y ~ z, data = correlation.df)


#Now that we have got the basics of R down, let's install phyloseq to be ready
#For Monday!

#Follow the instructions on the slides and give it a try!