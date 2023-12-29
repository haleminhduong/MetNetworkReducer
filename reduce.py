import tkinter as tk
from tkinter import filedialog, Listbox, Scrollbar, VERTICAL, messagebox
import os
import subprocess
import cobra
from cobra.io import read_sbml_model, load_matlab_model
import numpy as np
import pandas as pd
from numpy.linalg import matrix_rank


def load_model_file() -> None:
    global file_path
    selected_file_path: str = filedialog.askopenfilename(
            filetypes=[("MAT files", "*.mat"), ("SBML files", "*.xml"), ("All files", "*.*")])
    if selected_file_path:
        file_path = selected_file_path
        messagebox.showinfo("File Selected",
                            f"Model file selected: {file_path}")


def validate_user_input(delta_entry: tk.Entry,
                        bigM_entry: tk.Entry,
                        dof_entry: tk.Entry) -> list:
    try:
        default_delta, default_bigM, default_dof = 0.00001, 10000, 0
        delta: float = float(delta_entry.get()) if delta_entry.get()\
                else default_delta
        bigM: float = float(bigM_entry.get()) if bigM_entry.get()\
                else default_bigM
        dof: int = int(dof_entry.get()) if dof_entry.get() else default_dof
        return [delta, bigM, dof]
    except ValueError:
        messagebox.showerror("Input Error",
                             "Please enter valid numerical values for delta, \
                             bigm, and dof.")
        return None


def write_to_file(filename: str, data_list: list) -> None:
    with open(filename, 'w') as file:
        for item in data_list:
            file.write(f"{item}\n")


def process_input_file(input_file: str) -> list:
    with open(input_file, 'r') as file:
        lines = file.readlines()

    if ('no solution available\n' in lines):
        return ['infeasible']

    extracted_words: list = [word.replace('a$', '')
                             for line in lines
                             for word in line.split() 
                             if ('a$' in word)
                             and not ('<' in word)
                             and not ('v$' in word)
                             and not ('af$' in word)]
    return extracted_words


def write_to_prevSub(filename: str,
                     extracted_words: list,
                     iteration: int) -> None:
    with open(filename, 'a') as file:
        for word in extracted_words:
            file.write(f"{iteration} {word} 1\n")
        for word in [word for word in reactionsid
                     if word not in extracted_words]:
            file.write(f"{iteration} {word} 0\n")


def write_to_Iterations(filename: str, iteration: int) -> None:
    with open(filename, 'a') as file:
        file.write(f"{iteration}\n")


def check_DoF(minDoF: int,
              stoichiometric_matrix_df: pd.DataFrame,
              not_in_sub: list[str],
              num_sub_reactions: int) -> bool:
    matrix_pd_dropped = stoichiometric_matrix_df.drop(columns=not_in_sub)
    matrix_np = matrix_pd_dropped.to_numpy
    DoF: bool = num_sub_reactions - matrix_rank(matrix_np) >= minDoF
    return DoF


def run_pipeline() -> None:

    state: bool = False
    iteration: int = 0
    solutioncounter: int = 0

    with open(iterations_file, 'w') as file:
        file.write("0\n")

    while (not state):
        iteration += 1

        if os.path.exists(input_file):
            os.remove(input_file)

        # Update selections
        protected_reactions[:] = [reactionsid[i] for
                                  i in reactions_listbox.curselection()]
        protected_metabolites[:] = [metabolites[i] for
                                    i in metabolites_listbox.curselection()]

        params: list = validate_user_input(delta_entry, bigM_entry, dof_entry)
        if params is None:
            return

        # Writing to files
        write_to_file("zimpl_txt/Rev.txt", reversible_list)
        write_to_file("zimpl_txt/Irrev.txt", irreversible_list)
        write_to_file("zimpl_txt/l.txt", lower_bounds)
        write_to_file("zimpl_txt/u.txt", upper_bounds)
        write_to_file("zimpl_txt/PRxn.txt", protected_reactions)
        write_to_file("zimpl_txt/PMet.txt", protected_metabolites)
        write_to_file("zimpl_txt/params.txt", params)

        # Running subprocesses
        subprocess.run(['zimpl.exe', 'test.zpl'])
        subprocess.run(['scip.exe', '-f', 'test.lp', '-l', input_file])

        extracted_words: list[str] = process_input_file(input_file)


        if ('infeasible' in extracted_words):
            state = True
            return

        not_in_sub: list[str] = [word for word in reactionsid
                            if word not in extracted_words]

        num_sub_reactions: int = len(not_in_sub)

        minDoF: int = params[2]

        if (check_DoF(minDoF,
                      stoichiometric_matrix_df,
                      not_in_sub,
                      num_sub_reactions)):
            solutioncounter += 1
            write_to_prevSub(solution_file, extracted_words, solutioncounter)

        # Write output
        write_to_prevSub(output_file, extracted_words, iteration)
        write_to_Iterations(iterations_file, iteration)

    messagebox.showinfo("Pipeline Status",
                        "Pipeline executed with the provided parameters.")


def S_to_csv(stoichiometric_matrix: np.ndarray) -> pd.DataFrame:

    stoichiometric_matrix_df = pd.DataFrame(stoichiometric_matrix, index=internal_metabolites, columns=[rxn.id for rxn in model.reactions])

    # transform into readable format for ZIMPL
    stoichiometric_matrix_df_transposed = stoichiometric_matrix_df.T
    stoichiometric_matrix_df_reshaped = stoichiometric_matrix_df_transposed\
            .stack().reset_index()

    # output the matrix for ZIMPL
    stoichiometric_matrix_df_reshaped.to_csv("zimpl_txt/SMatrix_reshaped.txt",
                                         sep='\t')

    return stoichiometric_matrix_df


def add_functionality():
    global functionality_text_area
    global functionality_reactions_listbox
    # Function to update the text area with selected reactions
    selected_reactions = [reactionsid[i] for i in functionality_reactions_listbox.curselection()]
    formatted_reactions = ' + '.join([f"(1)v_{r}" for r in selected_reactions])
    functionality_text_area.insert(tk.END, formatted_reactions + " <= (1)\n")
    functionality_reactions_listbox.selection_clear(0, tk.END)


def get_functionality_text():
    global functionality_text_area, functionality_text
    # Get all text from the textarea
    functionality_text = functionality_text_area.get("1.0", tk.END)
    print(functionality_text)


def done_functionality():
    global functionality_window
    # Function to handle 'Done' button click
    # Extract and process data from functionality_text_area here
    get_functionality_text()
    functionality_window.withdraw()


# Function to open the new window for functionalities
def open_functionality_window():
    global functionality_window, functionality_text_area, functionality_reactions_listbox

    if functionality_window is not None:
        functionality_window.deiconify()

    else:

        functionality_window = tk.Toplevel(root)
        root.eval(f'tk::PlaceWindow {str(functionality_window)} center')

        # In the 'open_functionality_window' function
        info_button_functionality = tk.Button(functionality_window, text="Information", command=show_functionality_info)
        info_button_functionality.pack(pady=5)

        functionality_window.title("Functionalities")

        functionality_reactions_frame = tk.LabelFrame(functionality_window, text="Select Reactions for Functionality")
        functionality_reactions_frame.pack(fill="both", expand="yes", padx=10, pady=5)

        functionality_scrollbar = Scrollbar(functionality_reactions_frame, orient=VERTICAL)
        functionality_reactions_listbox = Listbox(functionality_reactions_frame, yscrollcommand=functionality_scrollbar.set, selectmode="multiple")
        functionality_scrollbar.config(command=functionality_reactions_listbox.yview)
        functionality_scrollbar.pack(side="right", fill="y")
        functionality_reactions_listbox.pack(side="left", fill="both", expand=True)

        for reaction in reactionsid:
            functionality_reactions_listbox.insert(tk.END, reaction)

        functionality_text_area = tk.Text(functionality_window, height=10, width=50)
        functionality_text_area.pack(padx=10, pady=10)

        add_functionality_button = tk.Button(functionality_window, text="Add new functionality", command=add_functionality)
        add_functionality_button.pack(pady=5)

        done_button = tk.Button(functionality_window, text="Done", command=done_functionality)
        done_button.pack(pady=5)

        functionality_window.protocol("WM_DELETE_WINDOW", done_functionality)


def show_main_info():
    messagebox.showinfo("Information", "Main Window Instructions: \n\n- Select protected reactions and metabolites.\n- Set parameters for delta, bigM, and DoF.\n- Add functionalities in 'Functionlities'\n- Click 'Run Pipeline' to execute.")


def show_functionality_info():
    messagebox.showinfo("Information", "Functionalities Instructions: \n\n- Select reactions to create a new functionality constraint and set the weights in the parentheses (1 by default).\n- Click 'Add new functionality' to start a new constraint.\n- Click 'Done' when all functionalities are added.")


# model = load_model("textbook")

# write_sbml_model(model, "test_e_coli.xml")
file_path: str = ''

load_model_file()

if ".mat" in file_path:
    model: cobra.Model = load_matlab_model(file_path)
else:
    # loads SBML network in
    model: cobra.Model = read_sbml_model(file_path)

# list for all the names of the metabolites
metabolites: list[str] = [metabolite.id for metabolite in model.metabolites]

# list for all the names of the reactions
reactionsid: list[str] = [reaction.id for reaction in model.reactions]

# create a stoichiometric matrix based on the network

num_metabolites: int = len(metabolites)
num_reactions: int = len(reactionsid)

#  Initialize an empty stoichiometric matrix
# stoichiometric_matrix: np.ndarray = np.zeros((num_metabolites, num_reactions))
#
#  Populate the stoichiometric matrix
# for i, reaction in enumerate(model.reactions):
#     for metabolite, coefficient in reaction.metabolites.items():
#         j = model.metabolites.index(metabolite)
#         stoichiometric_matrix[j, i] = coefficient

stoichiometric_matrix = cobra.util.array.create_stoichiometric_matrix(model)

# Identify external metabolites
external_metabolites = []
for met in model.metabolites:
    producing_reactions = [rxn for rxn in met.reactions if rxn.get_coefficient(met) < 0]
    consuming_reactions = [rxn for rxn in met.reactions if rxn.get_coefficient(met) > 0]

    if len(producing_reactions) == 0 or len(consuming_reactions) == 0:
        external_metabolites.append(met.id)

stoichiometric_matrix = cobra.util.array.create_stoichiometric_matrix(model)
internal_metabolites = [met.id for met in model.metabolites if met.id not in external_metabolites]
internal_met_index = [model.metabolites.index(met) for met in model.metabolites if met.id in internal_metabolites]
stoichiometric_matrix = stoichiometric_matrix[internal_met_index, :]

stoichiometric_matrix_df = S_to_csv(stoichiometric_matrix)

# list for all reactions
reactions: list[str] = [reaction.reaction for reaction in model.reactions]

# produce a list for all the reversible reactions
reversible_index: list[int] = ["=" in reaction for reaction in reactions]
reversible_list: list[str] = [reaction for reaction, reversible in
                   zip(reactionsid,
                       reversible_index) if reversible]

# produce a list for all the irreversible reactions
irreversible_list: list[str] = [irreversible for irreversible in reactionsid
                     if irreversible not in reversible_list]


# get lower and upper bounds
lower_bounds: list[float] = [reaction.lower_bound
                             for reaction in model.reactions]
upper_bounds: list[float] = [reaction.upper_bound
                             for reaction in model.reactions]


# set protected reactions and metabolites
protected_reactions: list[str] = []
protected_metabolites: list[str] = []


# set delta, bigM and DoF
params: list[float, float, int] = []

input_file: str = 'zimpl_txt/solution.txt'
output_file: str = 'zimpl_txt/prevMinimumSubnetworks.txt'
iterations_file: str = 'zimpl_txt/Iterations.txt'
solution_file: str = 'subnetworks.txt'
functionality_file: str = 'zimpl_txt/F.txt'
functionality_d_file: str = 'zimpl_txt/Fd.txt'

# Dmatrix: np.array = np.array([[0 for i in range(num_reactions)]])

Dmatrix: np.array = np.array([[0 for i in range(num_reactions)],
                              [0 for i in range(num_reactions)]])

biomassrxn_index = reactionsid.index("Growth")

Dmatrix[1][biomassrxn_index] = -1

Dmatrix = pd.DataFrame(Dmatrix)

Dmatrix.columns = [reaction.id for reaction in model.reactions]

Dmatrix_transposed = Dmatrix.T

Dmatrix_reshaped = Dmatrix_transposed.stack().reset_index()

Dmatrix_reshaped.to_csv("zimpl_txt/F.txt")

with open(solution_file, 'w') as file:
    pass

with open(output_file, 'w') as file:
    for rxn in reactionsid:
        file.write(f"0 {rxn} 0\n")

# Create the main window
root = tk.Tk()
root.eval('tk::PlaceWindow . center')
root.title("Metabolic Network Analysis Tool")
root.minsize(600, 400)

# In your main window setup
info_button_main = tk.Button(root, text="Information", command=show_main_info)
info_button_main.pack(pady=5)

# Listbox for selecting reactions
reactions_frame = tk.LabelFrame(root, text="Select Protected Reactions")
reactions_frame.pack(fill="both", expand="yes", padx=10, pady=5)
reactions_scrollbar = Scrollbar(reactions_frame, orient=VERTICAL)
reactions_listbox = Listbox(reactions_frame, yscrollcommand=reactions_scrollbar
                            .set, selectmode="multiple")
reactions_listbox.configure(exportselection=False)
reactions_scrollbar.config(command=reactions_listbox.yview)
reactions_scrollbar.pack(side="right", fill="y")
reactions_listbox.pack(side="left", fill="both", expand=True)
for reaction in reactionsid:
    reactions_listbox.insert(tk.END, reaction)

# Listbox for selecting metabolites
metabolites_frame = tk.LabelFrame(root, text="Select Protected Metabolites")
metabolites_frame.pack(fill="both", expand="yes", padx=15, pady=5)
metabolites_scrollbar = Scrollbar(metabolites_frame, orient=VERTICAL)
metabolites_listbox = Listbox(metabolites_frame,
                              yscrollcommand=metabolites_scrollbar.set,
                              selectmode="multiple")
metabolites_listbox.configure(exportselection=False)
metabolites_scrollbar.config(command=metabolites_listbox.yview)
metabolites_scrollbar.pack(side="right", fill="y")
metabolites_listbox.pack(side="left", fill="both", expand=True)
for metabolite in metabolites:
    metabolites_listbox.insert(tk.END, metabolite)

# Entries for delta and bigM
params_frame = tk.Frame(root)
params_frame.pack(fill="both", expand="yes", padx=15, pady=5)
tk.Label(params_frame, text="Delta:").pack(side="left")
delta_entry = tk.Entry(params_frame)
delta_entry.pack(side="left", padx=5)
tk.Label(params_frame, text="BigM:").pack(side="left")
bigM_entry = tk.Entry(params_frame)
bigM_entry.pack(side="left", padx=5)
tk.Label(params_frame, text="DoF:").pack(side="left")
dof_entry = tk.Entry(params_frame)
dof_entry.pack(side="left", padx=5)

# Button to open the functionality window
functionality_window = None
functionality_text = ""
functionalities_button = tk.Button(root, text="Functionalities", command=open_functionality_window)
functionalities_button.pack(pady=10)

# Button to run the pipeline
run_button = tk.Button(root, text="Run Pipeline", command=run_pipeline)
run_button.pack(pady=10)

# Start the GUI event loop
root.mainloop()
