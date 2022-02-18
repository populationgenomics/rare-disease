"""
class implementing bidirectional pedigree methods

1. Get the raw participant data, indexed by sample ID
2. Cast the pedigree members as objects,
    including immediate relationships
3. Assign children to members where appropriate
"""

from typing import Dict, List, Optional, Type
from dataclasses import dataclass
from csv import DictReader


@dataclass
class PedEntry:
    """
    class representing each participant from the PED
    """

    def __init__(
        self, fam, sample, father, mother, female, affected
    ):  # pylint: disable=too-many-arguments
        """

        :param fam:
        :param sample:
        :param father:
        :param mother:
        :param female:
        :param affected:
        """
        self.family: str = fam
        self.sample_id: str = sample

        # parental IDs as Strings or None
        self.father: str = father if father != '' else None
        self.mother: str = mother if mother != '' else None

        self.is_female = female == '2'
        self.affected = affected == '2'


@dataclass
class Participant:
    """
    dataclass representing a person within a family
    Type['Participant'] has to be used in order to have
    a self-referential class...

    Maybe that's a sign that it shouldn't be used...
    """

    details: PedEntry
    mother: Optional[Type['Participant']]
    father: Optional[Type['Participant']]
    children: List[Type['Participant']]


PED_KEYS = [
    '#Family ID',
    'Individual ID',
    'Paternal ID',
    'Maternal ID',
    'Sex',
    'Affected',
]


class PedigreeParser:
    """
    takes a PED file, and reads into a collection of linked-list like objects
    """

    def __init__(self, pedfile: str):
        """

        :param pedfile: path to a PED file
        """
        self.ped_dict = self.read_participants(pedfile)
        self.participants: Dict[str, Type['Participant']] = {}
        for sample_id in self.ped_dict.keys():
            self.populate_participants(sample_id=sample_id)
        self.apply_children()

    @staticmethod
    def read_participants(ped_file: str) -> Dict[str, PedEntry]:
        """
        Iterates through a PED file, and parses into a dict of Participants
        the dict is indexed by the Participant Sample ID string

        :param ped_file:
        :return:
        """

        with open(ped_file, 'r', encoding='utf-8') as handle:
            d_reader = DictReader(handle, delimiter="\t")
            participants = {
                party_line[PED_KEYS[1]]: PedEntry(
                    *[party_line.get(key) for key in PED_KEYS]
                )
                for party_line in d_reader
            }

        return participants

    def populate_participants(self, sample_id: str):
        """
        take a sample ID, and works backwards through the pedigree
        if we find a parent not already made into an object, create
        finally create a Participant for _this_ sample, including
        references to their parents, also as Participant objects
        :param sample_id:
        :return:
        """

        if sample_id in self.participants:
            return

        ped_sample = self.ped_dict.get(sample_id)
        if ped_sample.father is not None and ped_sample.father not in self.participants:
            self.populate_participants(ped_sample.father)
        if ped_sample.mother is not None and ped_sample.mother not in self.participants:
            self.populate_participants(ped_sample.mother)
        self.participants[sample_id] = Type['Participant'](
            details=ped_sample,
            mother=self.participants.get(ped_sample.mother),
            father=self.participants.get(ped_sample.father),
            children=[],
        )

    def apply_children(self):
        """
        iterates over the family members and adds children to nodes
        this allows for bi-directional inheritance checks
        :return:
        """

        # flick through all participants
        for participant in self.participants.values():
            # if this person has a father, add child to father
            if participant.father is not None:
                participant.father.children.append(participant)
            # repeat for mother
            if participant.mother is not None:
                participant.mother.children.append(participant)
