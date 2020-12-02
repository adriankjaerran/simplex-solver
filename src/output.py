""" Shows the output from running the simplex moethod
 



 - Inputted problem formulation
- Problem in standard form
    - If lower bounds - state substituted variable
    - If upper bound - state the variable and the three criterions that need to be checked during each iteration
    - Sentences describing which slack and artificial variables have been added and why (or if none is needed)
    - Whether or not a Phase 1 is needed

- Initial tableau
    - Initial basis variables

PER ITERATION
- Entering variable
- Exiting variable
- Direction

(finished)
- Optimal value
- Optimal values for variables
- Sensitivity analysis

INPUT FOR ADDING A NEW CONSTRAINT AND REOPTIMIZING
"""

# "The inputted problem formulation"

###################################
### Inputted problem formulation ###
###################################

#"Requirements for standard form"
# "- All constraints are equality constraints"
# "All variables must be non-negative (>=0)"


# "Step 1: Add slack variables and artificial variables to get the problem in standard form"

# "The following substituions were made:"



# "Lower bound substitutions"

#"Reference: p. 147-148"
#"Example: Assignment 5 Task 2"


#"x2: Since x2 has a lower bound of 5, it needs to be substituted to obtain the non-negativity structure
#" y2 = x2 - 5  -->  x2 = y2 + 5   , y2>=0"






#"Slack Variables"

#"Reference: p.81-82"
#"Example: " Example 4.1 p. 82

#"Constraint 1: A slack variable s1 is added to achieve equality constraint"
#" New constraint 1: 2x1 + 3x2 + s1 = 7 "





#"Artificial variables"

#"Reference: p.101-102"
#"Example:   p. 100 - 101"
#""          Assignment 5 Task 2" 

#"Constraint 2: An artificial variable a2 is added because the constraint is a >= or = constraint"
#"New constraint 2: 3x1 + x2 + a2 = 7"





#Free variables

#"Reference: p. 82"
#"Example:   Example 4.1 p. 82"

#"x3: is substituted by y3^+ - y3^- because x3 is a free variable 
#"x3 = y3+ - y3-    ,  y3+, y3- >= 0"





#STATE WHETHER PHASE 1 IS NEEDED 
#"Phase 1 is not needed because there are no artificial variables --> feasibility ensured"
#"Artificial variables introduced- Phase 1 is needed to indentify a initial feasible basic solution"



#Upper bound variables
#"Reference: p. 148"
#"Example:   Example 7.1  p.150"

#if no upper bounded variables
#   "No upper bounded variables, so the ordinary leaving variable criterion can be used"
#       "Leaving variable criterion: t^(k) = x_r^(k)/(-d_r^(k)) = min(over j) {x_j^(k)/(-d_j^(k)) | d_j^(k) < 0 }     (p.95)"

#if upper bounded variables
#   "Variable x2 is upper bounded by constraint 3. Therefore, the method simplex method for upper bounded variables is used"
    # Three criterions to consider:
#       1. t1: (Ordinary criterion) One of the current basic variables reaches its lower bound - this becomes leaving variable
#       2. t2: The entering variable reaches its upper bound - if so, it immediately becomes the leaving variable and the value is changed from 0 to the upper bound u_j
#       3. t3: One of the current basic variables reaches its upper bound - becomes leaving variable

#       Choose t^(k) = min{t1,t2,t3}
#       1. t1 = min(over j) {x_j^(k)/(-d_j^(k)) | d_j^(k) < 0}   (Ordinary criterion)     
#       2. t2 = u_p  (Entering variable reaches upper bound)
#       3. t3 = min {(u_j - x_j^(k)/d_j^(k)) | d_j^(k) > 0}   (The basic variable that first reaches its upper bound)
#   

# The resulting initial tableau

################################
### Problem in standard form ###
################################


# if(Phase 1 must be solved)
    # To find a feasible basis we have to solve the phase-1 problem
    # The goal is to minimize the sum of the artificial variables
    # min w = a1+a2 

    # If the optimal solution is w = 0, a feasible basic solution has been found --> ready for phase 2
    # if w > 0, the original problem has no feasible solution
    
    #####################
    ### Print tableau ###
    #####################

    #if upper bound method
        #The leaving basic variable is ..., since t^(iteration) = min{t1,t2,t3} = min {50,40,inf} = 40
        #where t1 = min{100/2, 80/1} = 50
        #      t2 = u1 = 40
        #      53 = min{upper bounds of basic variables} (often not defined --> inf:


    # Iteration 1, entering variable x1, leaving variable a2, distance vector, [1,0,-2,-3,0], #steps 5


    #####################
    ### Print tableau ###
    #####################

    # Iteration 2, entering variable x4, leaving variable a1, distance vector, [1,0,-2,-3,0], #steps 5

    # if(solution):
        # Found a feasible solution with x1=2, x2=0, x3=9, s1=9, s2=3
    
    # The objective function {objective_function} expressed in non_basic variables is
    # Objective function in non basic variables 
    
    ##############################################################
    ### Problem in standard form without artificial variables  ###
    ##############################################################


# Starting the phase-2 problem

# Iteration 1: Entering variable x1, leaving variable s2, distance vector [1,2,3,0,2], steps []
#               s2 is leaving because (reached lower bound / upper bound of leaving variable / the basis variable reached upper bound) 
# #if upper bound method
        #The leaving basic variable is ..., since t^(iteration) = min{t1,t2,t3} = min {50,40,inf} = 40
        #where t1 = min{100/2, 80/1} = 50
        #      t2 = u1 = 40
        #      53 = min{upper bounds of basic variables} (often not defined --> inf:
        
    #####################
    ### Print tableau ###
    #####################


# Iteration 2: Entering variable x2, leaving variable s1, distance vector [0,2,3,0,2], steps []
#               s2 is leaving because (reached lower bound / upper bound of leaving variable / the basis variable reached upper bound) 
#
    #####################
    ### Print tableau ###
    #####################



### Finished! The simplex method found an optimal solution ###

# Solution tableau 

    #####################
    ### Print tableau ###
    #####################

# Optimal solution: z^* = 37

# With variables with x1=2, x2=0, x3=9, s1=9, s2=3

# Sensitivity analysis
    # B^-1  
#
#



#INPUT: "Do you want to add a constraint and reoptimize using dual simplex? (Y/N)"
# or 
# INPUT "Do you want to change data and reoptimize using dual simplex? (Y/N)"

#If dual simplex
    #DUAL SIMPLEX

    #"Reference: p. 158-159"
    #"Example:   Example 7.3  p.160-162"

    # "Slack variable ... was added to the constraint --> equality --> normal form"
    # "The slack variable is added as a column to the initial tableau and is added to the basis"
    
    #Step 0: Start with the basic solution x^(0) which is dual feasible, i.e.

    #If Phase 1 is needed
        #Phase 1 is needed to find a initial dual feasible basic solution.
        #To be dual feasible the basic solution needs to satisfy
            #c_bar_j >= 0 (min problem)
            #c_bar_j <= 0 (max problem)
         
    #The new initial tableau is:

    #print new initial tableau

    #Step 1: Check the convergence criterion. The point x^(k) is an optimal solution if x_j^(k) >= 0 for all basic variables.
            #--> not fulfilled / Fulfilled

    #Step 2: Determine leaving basic variable. The new criterion is:
                # x_r^(k) = min (over j) {x_j^(k) | x_j^(k) < 0}
             #which gives that the basic variable x_r, represented on row s in the equation system, becomes the leaving basic variable.
             # -->Leaving basic variable = ...

    #Step 3: Determine entering basic variable. The new criterion is:
            # | C_bar_p / a_bar_sp | = min (over j) {|C_bar_j/a_bar_sj|  | A_bar_sj < 0}
            #which gives that x_p becomes the entering variable
            # --> Entering basic variable = ...
    
    #Step 4: Rewrite the system of equations as before, and compute x^(k+1).
            #Set k:= k+1 and go to step 1.
     


    #Iterations


    #(Finished)
    # Since all x_j >= 0 (no negative right hand side values) this solution is primal feasible. Since the solution to the dual and primal problem are equal, this solution is optimal.
    #Optimal solution:
        #z^* = ...
        #Variables = 
 

def printIntro(self):
    print("The inputted problem formulation")
    #Print the initial problem formulation
    self.printInputtedProblem()
    return self


def printInputtedProblem(self):
    print("###################################")
    print("### Inputted problem formulation ###")
    print("###################################")
    return self

def printSF(self):
    self.printSFIntro()
    self.printLowerBounds()
    self.printSlackVars()
    self.printArtificialVars()
    self.printFreeVars()
    self.printPhase1()
    self.printUpperBounds()
    self.printSFProblem()
    self.printTableauIntro()
    return self


def printSFIntro(self):
    print("Requirements for standard form:")
    print(" - All constraints are equality constraints")
    print(" - All variables must be non-negative (>=0) \n\n\n")
    print("The following substituions were made: \n\n")
    return self


def printLowerBounds(self):
    print("Lower bound substitutions \n\n")
    print("Reference: p. 147-148")
    print("Example: Assignment 5 Task 2 \n\n")

    
    #Template:
    #"x2: Since x2 has a lower bound of 5, it needs to be substituted to obtain the non-negativity structure
    #" y2 = x2 - 5  -->  x2 = y2 + 5   , y2>=0"
    
    return self

def printSlackVars(self):

    print("Slack Variables \n\n")

    print("Reference: p.81-82")
    print("Example:   Example 4.1 p. 82 \n\n")

    #Template
    #"Constraint 1: A slack variable s1 is added to achieve equality constraint"
    #" New constraint 1: 2x1 + 3x2 + s1 = 7 "

return self

def printArtificialVars(self):
    print("Artificial variables \n\n")

    print("Reference: p.101-102")
    print("Examples:   p. 100 - 101")
    print("           Assignment 5 Task 2")

    #Template
    #"Constraint 2: An artificial variable a2 is added because the constraint is a >= or = constraint"
    #"New constraint 2: 3x1 + x2 + a2 = 7"

    return self


def printFreeVars(self):


    print("Free variables \n\n")

    print("Reference: p. 82")
    print("Example:   Example 4.1 p. 82 \n\n")

    #Template
    #"x3: is substituted by y3^+ - y3^- because x3 is a free variable 
    #"x3 = y3+ - y3-    ,  y3+, y3- >= 0"



def printPhase1(self):
    if not self.require_two_phase:
        print ("Phase 1 is not needed because there are no artificial variables --> the initial basic solution is feasible")
    else :
        print("Artificial variables introduced --> Phase 1 is needed to indentify a initial feasible basic solution")

    return self       


def printUpperBounds(self):
    print("Upper Bounded variables \n\n")

    print("Reference: p. 148")
    print("Example:   Example 7.1  p.150 \n\n")

    if not self.upper_bounded:
        print("No upper bounded variables, so the ordinary leaving variable criterion can be used")
        print("     - Leaving variable criterion: t^(k) = x_r^(k)/(-d_r^(k)) = min(over j) {x_j^(k)/(-d_j^(k)) | d_j^(k) < 0 }     (p.95)"

    else:
    print("Variable ... is upper bounded by constraint ... . Therefore, the method simplex method for upper bounded variables is used \n\n")
    print("Three criterions to consider when choosing leaving variable:"
    print("     1. t1: (Ordinary criterion) One of the current basic variables reaches its lower bound - this becomes leaving variable")
    print("     2. t2: The entering variable reaches its upper bound - if so, it immediately becomes the leaving variable and the value is changed from 0 to the upper bound u_j")
    print("     3. t3: One of the current basic variables reaches its upper bound - becomes leaving variable \n\n")

    print("Choose t^(k) = min{t1,t2,t3}")
    print("     1. t1 = min(over j) {x_j^(k)/(-d_j^(k)) | d_j^(k) < 0}   (Ordinary criterion)")   
    print("     2. t2 = u_p                                              (Entering variable reaches upper bound)")
    print("     3. t3 = min {(u_j - x_j^(k)/d_j^(k)) | d_j^(k) > 0}      (The basic variable that first reaches its upper bound) \n\n")
    
    return self




    def printSFProblem(self):

        print("The problem in standard form \n\n")

        print("################################")
        print("### Problem in standard form ###")
        print("################################")

        return self


    def printTableauIntro(self):
        if self.require_two_phase:
            print(" \n\n PHASE 1")
            print("To find a feasible basis we have to solve the phase-1 problem. \n\n")
            print("The goal is to minimize the sum of the artificial variables.")
            print("     min w = a1+a2 \n\n")

            print("If the optimal solution is w = 0, a feasible basic solution has been found --> ready for phase 2")
            print("     if w > 0, the original problem has no feasible solution.")
            print(" \n\nThe resulting initial phase 1 tableau")   
        
        else:
            print("PHASE 2   (Ordinary simplex) \n\n")
            print(" \n\nThe resulting initial tableau")
        return self


    def printIteration(self, solved):
        self.printTableau()
        if self.require_two_phase:
            print("The leaving basic variable is ..., since t^(iteration) = min{t1,t2,t3} = min {... , ... , ...} = ...")
            print("t1 = min{100/2, 80/1} = 50")
            print("t2 = u1 = 40")
            print("t3 = min{upper bounds of basic variables} (often not defined --> inf)":

        if solved:
            self.printSolution()
            break

        self.printIterationValues()

        return self

    def printTableau(self):
        print("#####################")
        print("### Print tableau ###")
        print("#####################")
        print("\n\n")
        return self

    def printIterationValues(self):
        print("\n\n Iteration x \n\n")
        print("Entering variable: ...")
        print("Leaving variable: ...")
        print("Distance vector:  ...")

        return self

    def printSolution(self, dualSimplex):
        print("\n\n Finished! An optimal solution was found! \n\n")
        if self.require_two_phase:
            #Template
            print("\n\n Found a feasible solution with x1=2, x2=0, x3=9, s1=9, s2=3 \n\n")
            print("The objective function {objective_function} expressed in non_basic variables is")
            print("We now have the original problem in standard form")

            # kan vi bruke printSF? Hvordan er overgangen fra phase 1 til 2?

              ##############################################################
              ### Problem in standard form without artificial variables  ###
              ##############################################################
        
        if dualSimplex:
            
            print("Since all x_j >= 0 (no negative right hand side values) this solution is primal feasible. Since the solution to the dual and primal problem are equal, this solution is optimal.")
            print("Optimal solution:")
            print("     z^* = ...")
            print("     x1=... , x2 = ...") 



        else:    
            print("Optimal solution: z^* = ...")
            print("With variables ... = ... , ... = ... , ...=... ")
            
            #Template
            # With variables with x1=2, x2=0, x3=9, s1=9, s2=3)

            self.printSensitivityAnalysis()
            self.printContinueQuestion()

        return self

    def printSensitivityAnalysis(self):

        return self


    def printContinueQuestion(self):
        print("\n\nDo you want to add a constraint or a variable and reoptimize using dual simplex? (Y/N)")
        return self

    def printNextStepQuestion(self):
        print("Do you want to add a constraint (C) or add a variable (V)? (C/V)?")
        return self


    def printNextStepQuestion(self):
        print("Do you want to add a constraint (C) or add a variable (V)? (C/V)?")
        return self

 
    def printDSIntro(self, slack):
        print("\n\n DUAL SIMPLEX")
        print("Reference: p. 158-159")
        print("Example:   Example 7.3  p.160-162")

        if slack:
            print("Slack variable ... was added to the new constraint --> equality --> normal form")
            print("The slack variable extends the tableau with one column and one row as it is added to the basis")

    
        if self.require_two_phase:
            print("Phase 1 is needed to find a initial dual feasible basic solution.")
            print("To be dual feasible the basic solution needs to satisfy")
            print("     c_bar_j >= 0 (min problem)")
            print("     c_bar_j <= 0 (max problem)")

        
        print("Step 0: Start with the basic solution x^(0) which is dual feasible

        print("Step 1: Check the convergence criterion. The point x^(k) is an optimal solution if x_j^(k) >= 0 for all basic variables.")
                print("--> Fulfilled / Not fulfilled")

        print("\n\nStep 2: Determine leaving basic variable. The new criterion is:")
        print("     x_r^(k) = min (over j) {x_j^(k) | x_j^(k) < 0}")
        print("which gives that the basic variable x_r, represented on row s in the equation system, becomes the leaving basic variable.")
        print("     -->Leaving basic variable = ...")

        print("\n\nStep 3: Determine entering basic variable. The new criterion is:")
        print("     | C_bar_p / a_bar_sp | = min (over j) {|C_bar_j/a_bar_sj|  | A_bar_sj < 0}")
        print("which gives that x_p becomes the entering variable")
        print("     --> Entering basic variable = ...")
        
        print("\n\nStep 4: Rewrite the system of equations as before, and compute x^(k+1).")
        print("     Set k:= k+1 and go to step 1.")
        


        return self
