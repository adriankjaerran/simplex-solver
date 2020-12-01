import numpy as np


"""
Simplex solver

0) Parse input to numpy
1) Substitute lower bounds (xi_new = xi - l)
2) Add slack and artificial variables
3) Identity column for all basic variable columns

4.1) check if we need to do phase 1 simplex
    4.2 run simplex for phase 1 problem
    

5) Run simplex 

Simplex algorithm Pseudocode:
Input: list of basic variables, list of constraints, obj function, z-value
(self.basic_vars, A, b, c, z)

1) Find the largest reduced cost -> entering variable
    1.1) returns variable and direction for the largest reduced cost 
2) Check if optimal (no improving reduced cost) 
3) Find leaving variable: 
    step length t = min(t1, t2, t3)
    t1) Check number of steps before first basic variable reaches lower bound
    t2) t2 = upper bound of entering variable
    t3) Check number of steps before first basic variable reaches upper bound

4) New point: x_new = x + t*d
    Entering variable replaces leaving variable
    If t2: basis unchanged, entering = leaving
    if t2 or t3: substitute leaving variable xr_new = u - xr (u = upper bound)
    Perform pivot (identity column for all basic variable columns)
5) Go to step 1
6) Upon optimal solution: replace back substitutions
    xr = u - xr_new
    xi = xi_new + l
    




TODO implement the following methods:
    express_objective_in_nonbasic_variables(self,  )
    express_row_in_nonbasic_variabels(self, int row_num)
    pivot(self, entering_variable, leaving_variable)
    show_current_state(self, )
    to_standard_form(self, ) : should add the necessary slack and artificial variables
    





PRINT

- Inputted problem formulation
- Problem in standard form
- If lower bounds - state substituted variable




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


class LinearProgram:
    def __init__(self):

        self.upper_bounded = False  # if true, different check for leaving var
        self.phase_1 = False  # Whether phase 1 simplex must be performe

        # Original numpy arrays

        self.A_orig = None  # |vars| x |constraints|
        self.b_orig = None  # |constraints| x |1|
        self.c_orig = None  # |1| x |vars|
        self.z_orig = None  # |1| x |1|

        # Numpy arrays
        self.A = None  # |vars| x |constraints|
        self.b = None  # |constraints| x |1|
        self.c = None  # |1| x |vars|
        self.z = None  # |1| x |1|

        # Python lists/sets
        self.vars = None  # ['x1', 'x2', 'x3', 's1', 's2', 's3', 'a3']
        self.basic_vars = None  # ['x1', 'x2', 's1']    ordered by tableau
        self.nbasic_vars = None  # ['s2', 's3']
        self.slack_vars = None  # ['s1', 's2', 's3']
        self.artificial_vars = None  # ['a3']
        self.objective = None  # max or min
        self.free_vars = None  # ['x1', 'x2', 's1']

        # Transformations
        self.lb_substitutions = {}  # {'x2': ('y2', func, inv_func)}
        # self.fv_substitutions = {}  # {'x2': ('y2', func, inv_func)} # unsure if needed
        self.upper_bounds = {}  # {'x1': 9}

    def maximize(self, **obj_func):
        self.obj_func = obj_func
        self._objective = max
        return self

    def minimize(self, **obj_func):
        self.obj_func = obj_func
        self._objective = min
        return self

    def subject_to(self, *constraints, free_vars=None):
        self.constraints = list(constraints)
        self.free_vars = free_vars if free_vars else []
        return self

    def solve(self):
        # Load problem into numpy arrays
        self.c = np.array(list(self.obj_func.values()))

        # Load constraints into arrays

        # Convert to standard form
        self._to_standard_form()

        for const in self.constraints:
            print(const)
        print()
        print(self.A)
        print(self.b)
        # if self.phase_1:
        #   self._solve_phase_1()

        # Let's simplex
        pass

    def _to_standard_form(self):
        # TODO align changes with obj_func as well

        # Get all defined variables
        self.vars_set = set()
        for constraint in self.constraints:
            self.vars_set.update(constraint.keys())

        b = []
        remove_constraints = []

        # Create slack and artificial vars where needed
        for i, constraint in enumerate(self.constraints, start=1):
            # Extract data from constraints
            lhs = {k: v for k, v in constraint.items() if k not in ("op", "rhs")}
            op = constraint["op"]
            rhs = constraint["rhs"]

            slack_var = f"s{i}"
            artificial_var = f"a{i}"

            if op == "<=":
                # if only one variable: upper bounded = true
                if list(lhs.values()).count(0) == (len(lhs) - 1):
                    var = list(lhs.keys())[0]
                    self.upper_bounds[var] = rhs
                    remove_constraints.append(i)

                # add slack variable
                constraint[slack_var] = 1
                self.vars_set.add(slack_var)

            elif op == ">=":
                # if only one variable: lower bounded, perform substitution
                if len(lhs) == 1:
                    # replace the variable with y_i = x_i - rhs
                    var = list(lhs.keys())[0]
                    subst_var = f"y{i}"

                    # y_i = x_i - l --> x_i = y_i + l
                    self.lb_substitutions[var] = (
                        subst_var,
                        lambda x: x - rhs,
                        lambda y: y + rhs,
                    )
                    self.vars_set.discard(var)
                    constraint.pop(var)
                    self.vars_set.add(subst_var)
                    remove_constraints.append(i)

                # subtract slack variable, add artificial variable
                constraint[slack_var] = -1
                constraint[artificial_var] = 1
                self.vars_set.update([slack_var, artificial_var])

                # Flag phase 1 simplex
                self.phase_1 = True

            elif op == "=":
                # Add artificial variable
                constraint[artificial_var] = 1
                self.vars.add(artificial_var)
            b.append(rhs)

        for i in remove_constraints:
            print(self.constraints)
            self.constraints.pop(i - 1)
            b.pop(i - 1)

        # Handle free variables
        for fv in self.free_vars:
            # x_i = x_i+ - x_i-
            i = fv[-1]
            self.vars_set.discard(fv)

            for constraint in self.constraints:
                if fv in constraint.keys():
                    constraint[f"y{i}+"] = constraint[fv]
                    constraint[f"y{i}-"] = -1 * constraint[fv]
                    constraint.pop(fv)

        # Fill in default values if not present in constraint
        self.vars_set.update(set(self.obj_func.keys()))
        default_constraint = {v: 0 for v in self.vars_set}
        self.vars_set = set()
        for constraint in self.constraints:
            self.vars_set.update(constraint.keys())
            constraint = {**default_constraint, **constraint}

        # Create vector of upper bounds

        # Save constraints as np arrays
        A = []
        for constraint in self.constraints:
            lhs = {
                k: v for k, v in constraint.items() if k not in ("op", "rhs")
            }  # replace with .values() when smud order
            # TODO ensure smud order
            A.append(list(lhs.items()))
        self.A = np.array(A)
        self.b = np.array(b)

        return self

    def _single_iter(self):
        """
        Assumes the following:
            A: nxm matrix
        """
        pass

    #'x1'

    def _find_next_pivot(self):
        # Identify entering variable (most negative/positive entry)
        entering_col = np.argmax(-self.c if self.objective is max else self.c)

        # Criteria 1: Normal simplex criteria
        candidates = self.b / self.A[:, entering_col]
        candidates[candidates < 0] = float("inf")
        criteria_1 = np.min(candidates)

        # Criteria 2: Entering variable reaches upper bound
        criteria_2 = self.upper_bound.get(self.vars[entering_col], float("inf"))

        # Criteria 3: A current basic variable reaches its upper bounds
        candidates = self.b / self.A[:, entering_col]
        candidates[candidates < 0] = float("inf")
        criteria_1 = np.min(candidates)

        # Identify leaving variable

        # Identify direction vector
        pass

    def _pivot(self, leaving_var, entering_var):
        leaving_row = self.basic_vars.find(leaving_var)
        entering_col = self.vars.find(entering_var)

        # Update basis lists
        self.basic_vars[leaving_row] = entering_var

        # Normalize pivot
        pivot_value = self.A[leaving_row, entering_col]
        self.A[leaving_row, :] /= pivot_value

        # Row operations
        row_op_coeffs = self.A[:, entering_col]
        row_op_coeffs[leaving_row] = 0
        self.A -= self.A[leaving_row, :] * row_op_coeffs
        self.c -= self.A[leaving_row, :] * self.c[entering_col]

    ### Methods for extracting matrices and vectors ###

    ## Numpy arrays
    #   self.A = None                               # |vars| x |constraints|
    #   self.b = None                               # |constraints| x |1|
    #   self.c = None                               # |1| x |vars|
    #   self.z = None                               # |1| x |1|
    #
    #     #   # Python lists/sets
    #     #   self.vars = None                            # ['x1', 'x2', 'x3', 's1', 's2', 's3', 'a3']
    #     #   self.basic_vars = None                     # ['x1', 'x2', 's1']    ordered by tableau
    #   self.nbasic_vars = None                     # ['s2', 's3']
    #   self.slack_vars = None                      # ['s1', 's2', 's3']
    #   self.artificial_vars = None                 # ['a3']
    # self.objective = None                       # max or min

    # Matrices
    def B(self):
        # Find indexes of variables  currently in basis
        columns = [self.vars.index(i) for i in self.basic_vars]
        # Extract the indecies from the original coefficient matrix
        return self.A_orig[:, columns]

    def N(self):
        # Find indexes of variables  currently not in basis
        columns = [self.vars.index(i) for i in self.nbasic_vars]
        # Extract the indecies from the original coefficient matrix
        return self.A_orig[:, columns]

    def S(self):
        # Find indexes of variables  currently not in basis
        columns = [self.vars.index(i) for i in self.slack_vars]
        # Extract the indecies from the original coefficient matrix
        return self.A_orig[:, columns]

    # Coefficients
    def c(self):
        return self.c

    def c_B(self):
        # Find indecies of variables currently in basis
        columns = [self.vars.index(i) for i in self.basic_vars]
        return self.c[1, columns]

    def c_N(self):
        # Find indecies of variables currently not in basis
        columns = [self.vars.index(i) for i in self.nbasic_vars]
        return self.c[1, columns]

    def c_S(self):
        # Find indecies of variables in slack
        columns = [self.vars.index(i) for i in self.slack_vars]
        return self.c[1, columns]

    # Variables
    def x(self):
        return self.vars

    def x_B(self):
        return self.basic_vars

    def x_N(self):
        return self.nbasic_vars

    def x_S(self):
        return self.slack_vars

    # Values
    def b(self):
        return self.b.flatten()

    def z(self):
        return self.z.flatten()


_ = 0

## Problems
"""
lp.maximize(x1=3, x2=4).subject_to(
    ([1, _, 2],  '<=',   10),
    ([1, _, 2],  '<=',   10),
    ([1, _, 2],  '<=',   10),



"""
LinearProgram().maximize(x1=2, x2=3).subject_to(
    {"x1": 1, "x2": 1, "op": "<=", "rhs": 12},
    {"x1": 1, "x2": 2, "op": "<=", "rhs": 20},
    {"x1": 3, "x2": 1, "op": ">=", "rhs": 6},
    {"x1": 1, "x2": _, "op": ">=", "rhs": 1},
    {"x1": _, "x2": 9, "op": "<=", "rhs": 9},
).solve()
