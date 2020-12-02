import numpy as np
from copy import deepcopy
import utils

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
        self.z = np.zeros(None)  # |1| x |1|

        # Python lists/sets
        self.vars = []  # ['x1', 'x2', 'x3', 's1', 's2', 's3', 'a3']
        self.basic_vars = []  # ['x1', 'x2', 's1']    ordered by tableau
        self.objective = []  # max or min
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
        # Convert problem to mathematical formulation
        self._parse_problem()

        if self.require_two_phase:
            # Solve phase 1
            lp_phase_1 = self.copy()
            lp_phase_1.minimize(**{a: 1 for a in self.artificial_vars})
            #lp_phase_1.solve()

        # Let's simplex
        pass

    def _parse_problem(self):
        # Load problem into numpy arrays
        self.c = np.array(list(self.obj_func.values()))

        # Convert to standard form
        self._to_standard_form()
        pass

    def _to_standard_form(self):
        # TODO align changes with obj_func as well

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
                    self.require_two_phase = True

            elif op == "=":
                print(f"Adding {slack_var}")
                constraint[artificial_var] = 1
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
                    constraint["rhs"] -= constraint[key] * (constraint[key] - func(constraint[key])) # 
                    constraint.pop(key)

            # Substitute in objective function
            if key in self.obj_func.keys():
                self.obj_func[new_key] = self.obj_func[key] 
                constraint["rhs"] -= self.obj_func[key] * (self.obj_func[key] - func(self.obj_func[key])) # 
                self.obj_func.pop(key)


        # Handle free variables (x_i   -->   x_i+ - x_i-)
        for fv in self.free_vars:
            i = fv[-1]

            # Substitute in constraint
            for constraint in self.constraints:
                if fv in constraint.keys():
                    constraint[f"y{i}+"] = constraint[fv]
                    constraint[f"y{i}-"] = -1 * constraint[fv]
                    constraint.pop(fv)

            # Substitute in objective function
            if fv in self.obj_func.keys():
                self.obj_func[f"y{i}+"] = self.obj_func[fv]
                self.obj_func[f"y{i}-"] = -1 * self.obj_func[fv]
                self.obj_func.pop(fv)
                
        # Fill in default values if not present in constraint and objective function
        default_constraint = set(self.obj_func.keys())
        for constraint in self.constraints:
            default_constraint.update(set(constraint.keys()) - {"op", "rhs"})
        default_constraint = {v: 0 for v in default_constraint}

        # Add blanks in constraints and objective function
        self.obj_func = {**default_constraint, **self.obj_func}
        for i, constraint in enumerate(self.constraints):
            self.constraints[i] = {**default_constraint, **constraint}

        # Create variable set
        self.vars_set = set(default_constraint.keys())
        self.vars = list(self.vars_set)
        _x_vars = list(sorted(self.x_vars, key=lambda x: x[1:]))
        _s_vars = list(sorted(self.slack_vars))
        _a_vars = list(sorted(self.artificial_vars))
        self.vars = _x_vars + _s_vars + _a_vars

        # Create coefficient array
        _a_vars = list(sorted(self.artificial_vars))

        # Define objective function coefficients
        self.obj_func = {k: v for k, v in sorted(self.obj_func.items(), key=lambda tup: self.vars.index(tup[0]))}
        self.c = np.array(list(self.obj_func.values()))

        # Define matrices
        A = []
        b = []
        for constraint in self.constraints:
            # Sort constraint by order
            lhs = {k: v for k, v in constraint.items() if k not in ("op", "rhs")}  # replace with .values() when smud order1
            lhs = {k: v for k, v in sorted(lhs.items(), key=lambda tup: self.vars.index(tup[0]))}

            A.append(list(lhs.values()))
            b.append(constraint['rhs'])
        self.A = np.array(A)
        self.b = np.array(b)

        # Debugging
        print(f"vars = {self.vars}")
        print(f"free_vars = {self.free_vars}")
        print(f"nbasic_vars = {self.nbasic_vars}")
        print(f"slack_vars = {self.slack_vars}")
        print(f"artificial_vars = {self.artificial_vars}")
        
        print(f"\nc: {self.c}")
        print(f"\nbasic_vars = {self.basic_vars}")
        print(f"\nA: {self.A}")
        print(f"\nb: {self.b}")
        return self

    def _single_iter(self):
        """
        Assumes the following:
            A: nxm matrix
        """
        pass

    # 'x1'

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

    def copy(self):
        return deepcopy(self)

    # Matrices
    @property
    def B(self):
        b_indices = [self.vars.index(i) for i in self.basic_vars]
        return self.A_orig[:, b_indices]

    @property
    def N(self):
        nb_indices = [self.vars.index(i) for i in self.nbasic_vars]
        return self.A_orig[:, nb_indices]

    @property
    def S(self):
        s_indices = [self.vars.index(i) for i in self.slack_vars]
        return self.A_orig[:, s_indices]

    # Coefficients

    @property
    def c_B(self):
        b_indices = [self.vars.index(i) for i in self.basic_vars]
        return self.c[1, b_indices]

    @property
    def c_N(self):
        nb_indices = [self.vars.index(i) for i in self.nbasic_vars]
        return self.c[1, nb_indices]

    @property
    def c_S(self):
        s_indices = [self.vars.index(i) for i in self.slack_vars]
        return self.c[1, s_indices]

    # Variables
    def x(self):
        return self.vars

    def x_B(self):
        return self.basic_vars

    def x_N(self):
        return self.nbasic_vars

    def x_S(self):
        return self.slack_vars

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


_ = 0

# Problems
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
    # free_vars=['x2']
).solve()
