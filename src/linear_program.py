import numpy as np
import fractions
np.set_printoptions(edgeitems=30, linewidth=1000)
np.set_printoptions(formatter={'float':lambda x: f"{str(fractions.Fraction(x).limit_denominator()) if np.isfinite(x) else '∞':4}"})


from copy import deepcopy
import utils

GLOBAL_COUNTER = 0

class LinearProgram:
    def __init__(self):

        self.upper_bounded = False  # if true, different check for leaving var
        self.require_two_phase = False  # Whether phase 1 simplex must be performe

        # Original numpy arrays
        self.A_orig = np.zeros(None)  # |vars| x |constraints|
        self.b_orig = np.zeros(None)  # |constraints| x |1|
        self.c_orig = np.zeros(None)  # |1| x |vars|
        self.z_orig = np.zeros(None)  # |1| x |1|

        # Numpy arrays
        self.A = np.zeros(None)  # |vars| x |constraints|
        self.b = np.zeros(None)  # |constraints| x |1|
        self.c = np.zeros(None)  # |1| x |vars|
        self.z = np.zeros((1, 1))  # |1| x |1|

        # Python lists/sets
        self.vars = []  # ['x1', 'x2', 'x3', 's1', 's2', 's3', 'a3']
        self.basic_vars = []  # ['x1', 'x2', 's1']    ordered by tableau
        self.objective = None  # max or min
        self.obj_func = {}  # {'x1': 2, 'x2': 3}
        self.free_vars = []  # ['x1', 'x2', 's1']

        # AUTO-GENERATED PROPERTIES
        # self.nbasic_vars = None  # ['s2', 's3']
        # self.slack_vars = None  # ['s1', 's2', 's3']
        # self.artificial_vars = None  # ['a3']

        # Transformations
        self.lb_substitutions = {}  # {'x2': ('y2', func, inv_func)}
        # self.fv_substitutions = {}  # {'x2': ('y2', func, inv_func)} # unsure if needed
        self.upper_bounds = {}  # {'x1': 9}

    def maximize(self, **obj_func):
        self.obj_func = obj_func
        self.objective = max
        return self

    def minimize(self, **obj_func):
        self.obj_func = obj_func
        self.objective = min
        return self

    def subject_to(self, *constraints, free_vars=None):
        self.constraints = list(constraints)
        self.free_vars = free_vars if free_vars else []
        return self

    def solve(self):
        if self.require_two_phase:
            # Solve phase 1
            lp_phase_1 = self.copy()
            lp_phase_1.require_two_phase = False
            lp_phase_1.minimize(**{a: 1 for a in lp_phase_1.artificial_vars}) #minimize(x1=2)

            # Formulate objective function
            lp_phase_1.obj_func = {**lp_phase_1.default_constraint, **lp_phase_1.obj_func}
            lp_phase_1.obj_func = utils.reorder_by_index(lp_phase_1.obj_func, lp_phase_1.vars)
            lp_phase_1.c = -np.array(list(lp_phase_1.obj_func.values()), dtype="float64").reshape(1, -1)
            lp_phase_1.z = np.zeros((1, 1))
            
            pairs = list(enumerate(lp_phase_1.basic_vars))
            artificial_rows = np.take(lp_phase_1.A, [i for i, v in pairs if v.startswith('a')], axis=0)
            artificial_rhs = np.take(lp_phase_1.b, [i for i, v in pairs if v.startswith('a')], axis=0)
            lp_phase_1.c += np.sum(artificial_rows, axis=0)
            lp_phase_1.z += np.sum(artificial_rhs, axis=0)
            
            print("\n"*5, "LP PHASE 1!")
            print(lp_phase_1)
            lp_phase_1.solve()

            # Ensure that phase 1 is solvable
            if lp_phase_1.z > 0 :
                print("infeasible!!!")
                return

            # Update variables of current instance
            self.A = lp_phase_1.A
            self.b = lp_phase_1.b
            self.basic_vars = lp_phase_1.basic_vars

            # Remove artificial variables
            self.A = self.A[:, :-len(self.artificial_vars)]
            self.c = self.c[:, :-len(self.artificial_vars)]
            self.vars = self.x_vars + self.slack_vars

            # express z with non-basic vars
            self.z = self.z - np.sum(self.b * self.c_B.reshape(-1, 1), axis=0) 
            self.c = self.c - np.sum(self.A * self.c_B.reshape(-1, 1), axis=0) 
            # print(100*"#")
            # print(self)

        # Let's simplex
        for i in range(10):
            print("\n"*5, f"ITERATION {i}\n")
            print(self)
            self._single_iter()
            print("\n"*15)

            
            # Check for optimality
            if np.all((self.c >= 0) if self.objective is max else (self.c <= 0)):
                print("\nSOLUTION FOUND!!")
                print(self)
                return
            # entering_idx = np.argmax(-self.c if self.objective is max else self.c)
            # entering_var = self.vars[entering_idx]

    def compile(self):
        """Converts the problem to mathematical formulation."""
        # Get all defined variables
        self.vars_set = set()
        for constraint in self.constraints:
            self.vars_set.update(constraint.keys())

        b = []
        remove_constraints = set()

        # Create slack and artificial vars where needed
        for i, constraint in enumerate(self.constraints):
            # Extract data from constraints
            lhs = {k: v for k, v in constraint.items() if k not in ("op", "rhs")}
            op = constraint["op"]
            rhs = constraint["rhs"]

            slack_var = f"s{i+1}"
            artificial_var = f"a{i+1}"

            if op == "<=":
                # if only one variable: upper bounded = true
                if utils.is_variable_constraint(lhs):
                    upper_bound = utils.get_var_from_constraint(lhs)
                    self.upper_bounds[upper_bound] = rhs
                    remove_constraints.add(i)
                else:
                    print(f"Adding {slack_var}")
                    constraint[slack_var] = 1
                    self.basic_vars.append(slack_var)

            elif op == ">=":
                if utils.is_variable_constraint(lhs):
                    # Replace the variable with y_i = x_i - rhs
                    var = utils.get_var_from_constraint(lhs)
                    print("VAR", var)
                    print("RHS", rhs)
                    self.lb_substitutions[var] = (
                        f"y{var[1:]}",
                        (lambda rhs: lambda x: x - rhs)(rhs),
                        (lambda rhs: lambda y: y + rhs)(rhs),
                    )
                    constraint.pop(var)
                    remove_constraints.add(i)
                else:
                    print(f"Adding {slack_var}, {artificial_var}")
                    constraint[slack_var] = -1
                    constraint[artificial_var] = 1
                    self.basic_vars.append(artificial_var)
                    self.require_two_phase = True

            elif op == "=":
                print(f"Adding {artificial_var}")
                constraint[artificial_var] = 1
                self.basic_vars.append(artificial_var)

            b.append(rhs)

            # Update actual constraint list with modified constraint
            self.constraints[i] = constraint

        # Remove substituted singular >= - constraints
        self.constraints = [con for i, con in enumerate(self.constraints) if i not in remove_constraints]
        b = [value for i, value in enumerate(b) if i not in remove_constraints]

        # Remove substituted variables


        # Execute substitutions
        for key, values in self.lb_substitutions.items():
            new_key, func, _ = values
            print(f"Substituting {key} with {new_key}")

            # Substitute in constraint
            for constraint in self.constraints:
                if key in constraint.keys():
                    constraint[new_key] = constraint[key] 
                    constraint["rhs"] += constraint[key] * func(0)
                    constraint.pop(key)

            # Substitute in objective function
            if key in self.obj_func.keys():
                self.obj_func[new_key] = self.obj_func[key] 
                self.z += self.obj_func[key] * -func(0)
                self.obj_func.pop(key)

        # Handle free variables (x_i   -->   x_i+ - x_i-)
        for fv in self.free_vars:
            i = fv[-1]

            # Substitute in constraint
            for constraint in self.constraints:
                if fv in constraint.keys():
                    new_plus = f"y{i}+"
                    new_minus = f"y{i}-"
                    print(f"Substituting {fv} with {new_plus} and {new_minus}")

                    constraint[new_plus] = constraint[fv]
                    constraint[new_minus] = -1 * constraint[fv]
                    constraint.pop(fv)

            # Substitute in objective function
            if fv in self.obj_func.keys():
                self.obj_func[f"y{i}+"] = self.obj_func[fv]
                self.obj_func[f"y{i}-"] = -1 * self.obj_func[fv]
                self.obj_func.pop(fv)
                
        # Fill in default values if not present in constraint and objective function
        self.default_constraint = set(self.obj_func.keys())
        for constraint in self.constraints:
            self.default_constraint.update(set(constraint.keys()) - {"op", "rhs"})
        self.default_constraint = {v: 0 for v in self.default_constraint}

        # Add blanks in constraints and objective function
        self.obj_func = {**self.default_constraint, **self.obj_func}
        for i, constraint in enumerate(self.constraints):
            self.constraints[i] = {**self.default_constraint, **constraint}

        # Create variable set
        self.vars_set = set(self.default_constraint.keys())
        self.vars = list(self.vars_set)
        _x_vars = list(sorted(self.x_vars, key=lambda x: x[1:]))
        _s_vars = list(sorted(self.slack_vars))
        _a_vars = list(sorted(self.artificial_vars))
        self.vars = _x_vars + _s_vars + _a_vars

        # Create upper bound vector
        default_bound = {k: float('inf') for k in self.vars_set}
        self.upper_bounds = {**default_bound, **self.upper_bounds}
        self.upper_bounds = utils.reorder_by_index(self.upper_bounds, self.vars)
        self.upper_bounds_vec = np.array(list(self.upper_bounds.values())).reshape(1, -1)

        # Create coefficient array
        _a_vars = list(sorted(self.artificial_vars))

        # Define objective function coefficients
        self.obj_func = utils.reorder_by_index(self.obj_func, self.vars)
        self.c = -np.array(list(self.obj_func.values()), dtype="float64").reshape(1, -1)

        # Define matrices
        A, b = [], []
        for constraint in self.constraints:
            # Sort constraint by order
            lhs = {k: v for k, v in constraint.items() if k not in ("op", "rhs")}
            lhs = utils.reorder_by_index(lhs, self.vars)
            A.append(list(lhs.values()))
            b.append(constraint['rhs'])
        self.A = np.array(A, dtype="float64")
        self.b = np.array(b, dtype="float64").reshape(-1, 1)
        return self

    def _single_iter(self):
        criteria, leaving_var, entering_var = self._find_next_pivot()
        self._pivot(criteria, leaving_var, entering_var)

    def _find_next_pivot(self):
        # Identify entering variable (most negative/positive entry)
        #entering_idx = np.argmin(-self.c if self.objective is max else self.c)
        entering_idx = np.argmin(self.c if self.objective is max else -self.c)
        entering_var = self.vars[entering_idx]

        #TODO: output lists we are minimzing

        # Criteria 1: Normal simplex criteria
        numerator = self.b.ravel()
        denominator = self.A[:, entering_idx].ravel()
        fraction = np.divide(numerator, denominator, out=np.ones_like(denominator)*float('inf'), where=denominator>0)
        criteria_1 = fraction.min()
        criteria_1_var = self.basic_vars[np.argmin(fraction)]

        # Criteria 2: Entering variable reaches upper bound
        criteria_2 = self.upper_bounds_vec[0, entering_idx]
        criteria_2_var = self.vars[entering_idx]

        # Criteria 3: A current basic variable reaches its upper bounds
        numerator = np.take(self.upper_bounds_vec, self._basic_indices) - self.b.ravel()
        denominator = self.A[:, entering_idx]
        fraction = np.divide(numerator, -denominator, out=np.ones_like(denominator)*float('inf'), where=denominator<0)
        criteria_3 = fraction.min() 
        criteria_3_var = self.basic_vars[np.argmin(fraction)]
        
        # Identify leaving variable
        defining_criteria = np.argmin([criteria_1, criteria_2, criteria_3])
        leaving_var = [criteria_1_var, criteria_2_var, criteria_3_var][defining_criteria]

        # Identify direction vector
        direction = np.zeros(len(self.vars))
        direction[entering_idx] = 1
        for i, var in enumerate(self.basic_vars):
            direction[self.vars.index(var)] = -self.A[i, entering_idx]
        
        print(f"Direction: {tuple(self.vars)} = {direction}")
        print(f"T-values: {np.array([criteria_1, criteria_2, criteria_3])}")
        print(f"T-variables: {np.array([criteria_1_var, criteria_2_var, criteria_3_var])}")
        print(f"Defining criteria: {defining_criteria} \nEntering var: {entering_var} \nLeaving var: {leaving_var}")

        return defining_criteria, leaving_var, entering_var

    def _pivot(self, criteria, leaving_var, entering_var):
        # Check if we need to substitute
        if criteria == 1:
            col_idx = self.vars.index(leaving_var)

            # Update problem variables
            if not leaving_var.endswith("^"):
                sub_var = f"{leaving_var[:2]}^"
            else: 
                sub_var = f"{leaving_var[:2]}"
            self.vars[col_idx] = sub_var
            print(f"Substituting {leaving_var} with {sub_var}")

            self.A[:, col_idx] *= -1
            self.c[:, col_idx] *= -1

            self.b += self.upper_bounds_vec.ravel()[col_idx] * self.A[:, col_idx].reshape(-1, 1)
            self.z += self.upper_bounds_vec.ravel()[col_idx] * self.c[:, col_idx]

            # print(f'#### Result after iteration ')
            # print(f"\nbasic_vars = {self.basic_vars}")
            # print(f"\nc: {self.c}")
            # print(f"\nA: {self.A}")
            # print(f"\nb: {self.b}")
            # print(f"\nz: {self.z}\n\n\n")
            return 

        elif criteria == 2:
            row_idx = self.basic_vars.index(leaving_var)
            col_idx = self.vars.index(leaving_var)

            # Update problem variables
            leaving_var = f"{leaving_var[:2]}^"
            self.vars[col_idx] = leaving_var
            self.basic_vars[row_idx] = leaving_var

            # Row operations
            self.A[row_idx, :] *= -1
            self.A[row_idx, col_idx] = 1

            # Modify RHS
            self.b[row_idx] = -self.b[row_idx] + self.upper_bounds_vec[0, col_idx]

        # Perform pivot
        # Update basis lists
        row_idx = self.basic_vars.index(leaving_var)
        self.basic_vars[row_idx] = entering_var

        # Normalize pivot
        col_idx = self.vars.index(entering_var)
        pivot_value = self.A[row_idx, col_idx]
        self.A[row_idx, :] = self.A[row_idx, :] / pivot_value
        self.b[row_idx, :] = self.b[row_idx, :] / pivot_value

        # Row operations
        row_op_coeffs = self.A[:, col_idx].copy()
        row_op_coeffs[row_idx] = 0
    
        self.A -= self.A[row_idx, :] * row_op_coeffs.reshape(-1, 1)
        self.b -= self.b[row_idx] * row_op_coeffs.reshape(-1, 1)

        # Objective function operations
        self.z -= self.b[row_idx] * self.c[0, col_idx]
        self.c -= self.A[row_idx, :].reshape(1, -1) * self.c[0, col_idx]
        
        # print(f"#### Result after iteration ####")
        # print(f"\nbasic_vars = {self.basic_vars}")
        # print(f"\nc: {self.c}")
        # print(f"\nA: {self.A}")
        # print(f"\nb: {self.b}")
        # print(f"\nz: {self.z}")
        # print(self)

    def copy(self):
        return deepcopy(self)

    def __str__(self):
         # Debugging
        output = ""
        output += f"vars = {self.vars}\n"
        output += f"free_vars = {self.free_vars}\n"
        output += f"nbasic_vars = {self.nbasic_vars}\n"
        output += f"slack_vars = {self.slack_vars}\n"
        output += f"artificial_vars = {self.artificial_vars}\n"
        output += f"basic_vars = {self.basic_vars}\n"
        output += f"upper_bounds_vec = {self.upper_bounds_vec}\n"
        # output += f"c: {self.c}\n"
        # output += f"A: {self.A}\n"
        # output += f"b: {self.b}\n"
        # output += f"z: {self.z}\n"
        output += str(np.hstack((np.vstack((np.ones((1, 1)), np.zeros_like(self.b))), np.vstack((self.c, self.A)), np.vstack((self.z, self.b)))))

        return output

    # Matrices
    @property
    def B(self):
        return self.A_orig[:, self._basic_indices]

    @property
    def N(self):
        return self.A_orig[:, self._nbasic_indices]

    @property
    def S(self):
        return self.A_orig[:, self._slack_indices]

    # Coefficients

    @property
    def c_B(self):
        return self.c[0, self._basic_indices]

    @property
    def c_N(self):
        return self.c[0, self._nbasic_indices]

    @property
    def c_S(self):
        return self.c[0, self._slack_indices]

    # Variables
    def x(self):
        return self.vars

    def x_B(self):
        return self.basic_vars

    def x_N(self):
        return self.nbasic_vars

    def x_S(self):
        return self.slack_vars

    # Indices for variable sets
    @property
    def _basic_indices(self):
        return np.array([self.vars.index(i) for i in self.basic_vars])
        
    @property
    def _nbasic_indices(self):
        return np.array([self.vars.index(i) for i in self.nbasic_vars])
        
    @property
    def _slack_indices(self):
        return np.array([self.vars.index(i) for i in self.slack_vars])
        
    @property
    def _artificial_indices(self):
        return np.array([self.vars.index(i) for i in self.artificial_vars])

    # Properties
    @property
    def x_vars(self):
        return [v for v in self.vars if v.startswith('x') or v.startswith('y')]

    @property
    def nbasic_vars(self):
        return [v for v in self.vars if v not in self.basic_vars]

    @property
    def slack_vars(self):
        return [v for v in self.vars if v.startswith('s')]

    @property
    def artificial_vars(self):
        return [v for v in self.vars if v.startswith('a')]




# Problems
"""
lp.maximize(x1=3, x2=4).subject_to(
    ([1, _, 2],  '<=',   10),
    ([1, _, 2],  '<=',   10),
    ([1, _, 2],  '<=',   10),



"""
_ = 0
x1 = "x1"
x2 = "x2"
x3 = "x3"
x4 = "x4"
x5 = "x5"
op = "op"
rhs = "rhs"

# Assignment 5
lp1 = LinearProgram().maximize(x1=3, x2=4).subject_to(
    {x1: 1, x2: 1, op: "<=", rhs: 12},
    {x1: 1, x2: 2, op: "<=", rhs: 20},
    {x1: 3, x2: 1, op: ">=", rhs: 6},
    {x1: 1, x2: _, op: ">=", rhs: 1},
    {x1: _, x2: 9, op: "<=", rhs: 9},
    # free_vars=['x2']
)

# LRV p. 150
lp2 = LinearProgram().maximize(x1=30, x2=20).subject_to(
    {x1: 2, x2: 1, op: "<=", rhs: 100},
    {x1: 1, x2: 1, op: "<=", rhs: 80},
    {x1: 1, x2: _, op: "<=", rhs: 40},
    {x1: _, x2: 1, op: "<=", rhs: 30},
)

#LRV p.96
lp3 =  LinearProgram().maximize(x1=30, x2=20).subject_to(
    {x1: 2, x2: 1, op: "<=", rhs: 100},
    {x1: 1, x2: 1, op: "<=", rhs: 80},
    {x1: 1, x2: _, op: "<=", rhs: 40},
)

# Exam Fall 2018 exercise 2
lp4 = LinearProgram().maximize(x1=6, x2=-4, x3=1).subject_to(
    {x1: 1, x2: -2, x3: _, op: "<=", rhs: 0},
    {x1: 2, x2: 1, x3: 1, op: "<=", rhs: 10},
    {x1: 1, x2: _, x3: _, op: "<=", rhs: 3},
    {x1: _, x2: 1, x3: _, op: "<=", rhs: 2},
    {x1: _, x2: _, x3: 1, op: "<=", rhs: 5},
)

# Exam Fall 2019 exercise 1
lp5 =  LinearProgram().maximize(x1=3, x2=-1).subject_to(
    {x1: 1, x2: 1, op: ">=", rhs: 2},
    {x1: 1, x2: 4, op: ">=", rhs: 4},
    {x1: -1, x2: 2, op: "<=", rhs: 4},
    {x1: 3, x2: 1, op: "<=", rhs: 6},
)

# Assignment 3
lp6 = LinearProgram().maximize(x1=2, x2=1).subject_to(
    {x1: 1, x2: -3, op: "<=", rhs: 3},
    {x1: -1, x2: 1, op: "<=", rhs: 4},
    {x1: 1, x2: 2, op: "<=", rhs: 8},
    # free_vars=['x2']
)


# Assignment 10, exercise 1

lp7 = LinearProgram().maximize(x1=3, x2=2).subject_to(
    {x1: 1, x2: 2, op: ">=", rhs: 2},
    {x1: 1, x2: -2, op: "<=", rhs: 1},
    {x1: 1, x2: 1, op: "<=", rhs: 5}
)

# Assignment 11, exercise 2

lp8 = LinearProgram().maximize(x1=4, x2=8, x3 = 9, x4=10).subject_to(
    {x1: 3, x2: 5, x3:6, x4: 8, op: "<=", rhs: 10},
    {x1: 1, op:"<=", rhs: 1},
    {x2: 1, x3: 1, x4: 1, op:"<=", rhs: 1}
)

lp9 = LinearProgram().maximize(x1=4, x2=8, x4=10).subject_to(
    {x1: 3, x2: 5, x4: 8, op: "<=", rhs: 10},
    {x1: 1, op:"<=", rhs: 1},
    {x2: 1, op:"<=", rhs: 0},
    #{x3: 1, op:"<=", rhs: 1},
    #{x4: 1, op:"<=", rhs: 1},
)


lp8.compile().solve()