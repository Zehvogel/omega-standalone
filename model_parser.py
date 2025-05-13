import re
from math import sqrt, pi


parameter_pattern = re.compile(r"parameter\s+(\w+)\s+=\s+(\d+\.?\d*(?:[eE]-?\d+)?)")
derived_pattern = re.compile(r"derived\s+(\w+)\s+=\s+([\w\s\(\)\-\+\*\/\.]*[\w\)])")


def get_parsed_parameters(path: str) -> dict[str, float]:
    results: dict[str, float] = {}
    with open(path) as file:
        text = file.read()
    for match in re.finditer(parameter_pattern, text):
        name, value = match.group(1, 2)
        # print(name, float(value))
        results[name] = float(value)
    for match in re.finditer(derived_pattern, text):
        name, value = match.group(1, 2)
        # don't run this on random model files pls :)
        derived_value = eval(value, (results | {"sqrt": sqrt, "pi": pi}))
        print(name, derived_value)
        results[name] = derived_value

    return results


if __name__  == "__main__":
    get_parsed_parameters("SM_ac.mdl")