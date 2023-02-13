import xml.etree.ElementTree as ET
import argparse

class XMLParser:
    def __init__(self, xml_file):
        self.xml_file = xml_file

    def parse(self):
        tree = ET.parse(self.xml_file)
        root = tree.getroot()

        result = {}
        for alignment in root.findall('./alignment'):
            for sequence in alignment.findall('./sequence'):
                seq_name = sequence.find('./seq-name').text
                seq_data = sequence.find('./seq-data').text
                result[seq_name] = seq_data

        return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract dictionary from XML file')
    parser.add_argument('-x', '--xml_file', help='Path to the XML file')
    parser.add_argument('-f', '--fasta', help='Print the output in FASTA format')
    args = parser.parse_args()

    xml_parser = XMLParser(args.xml_file)
    result = xml_parser.parse()

    if args.fasta:
        with open(args.fasta, 'w') as fi:
            for key, value in result.items():
                fi.write(f'>{key}\n{value}')
    else:
        print(result)
