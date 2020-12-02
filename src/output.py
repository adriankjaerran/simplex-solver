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
"""

# The inputed problem formulation

###################################
### Inputed problem formulation ###
###################################

#Requirements for standard form
# "- All constraints are equality constraints"
# "All variables must be non-negative (>=0)"


# Step 1: Add slack variables and artificial variables to get the problem in standard form

# The following substituions were made



# Lower bound substitutions 

#"Reference: p. 147-148"
#"Example: Assignment 5 Task 2"


#"x2: Since x2 has a lower bound of 5, it needs to be substituted to obtain the non-negativity structure
#" y2 = x2 - 5  -->  x2 = y2 + 5   , y2>=0"






# Slack Variables

#"Reference: p.81-82"
#"Example: " Example 4.1 p. 82

#"Constraint 1: A slack variable s1 is added to achieve equality constraint"
#" New constraint 1: 2x1 + 3x2 + s1 = 7 "





#Artificial variables

#"Reference: p.101-102"
#"Example:   p. 100 - 101"
#            Assignment 5 Task 2   " 

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
        #      53 = min{upper bounds of basic variables} (often not defined --> inf)


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
        #      53 = min{upper bounds of basic variables} (often not defined --> inf)
        
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





 



basis = ["x2","x3"