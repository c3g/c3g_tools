#!/usr/bin/env python3

"""
Merges sub-yaml files within template yaml provided.
"""

# General import
import argparse
import sys
import yaml

def parseoptions():
    """Command line options"""
    parser = argparse.ArgumentParser(description="Merges sub-yaml files within template yaml provided.")

    parser.add_argument('-t',
                        '--template',
                        help="Template yaml file",
                        required=True)
    parser.add_argument('-y',
                        '--sub_yaml',
                        help="Section yaml file(s) to be merged within the template.",
                        required=True,
                        nargs='+')
    parser.add_argument('-o',
                        '--output',
                        help="Output yaml file name",
                        required=False)

    return parser.parse_args()


def main():
    """main function"""

    # ARGS
    args = parseoptions()

    with open(args.template, "r", encoding="utf8") as template_file:
        template = yaml.safe_load(template_file)

    for file_name in args.sub_yaml:
        with open(file_name, 'r', encoding="utf8") as yaml_file:
            sub_yaml = yaml.safe_load(yaml_file)
        for index, item in enumerate(template["module_order"]):
            try:
                if item in sub_yaml.keys():
                    template["module_order"][index] = sub_yaml[item]
            except:
                pass

    if args.output:
        with open(args.output, 'w') as output_file:
            yaml.dump(template, output_file)
    else:
        yaml.dump(template, sys.stdout)

if __name__ == "__main__":
    main()
