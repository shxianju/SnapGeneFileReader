import struct
import xmltodict
import os
import time
import textwrap

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
import html2text
import re

import json


html_parser = html2text.HTML2Text()
html_parser.ignore_emphasis = True
html_parser.ignore_links = True
html_parser.body_width = 0
html_parser.single_line_break = True


def wrap_text(text, leading_space):
    wrapper = textwrap.TextWrapper()
    wrapper.width = 80
    wrapper.subsequent_indent = " " * leading_space
    text = text.rstrip("\n")
    text_lines = text.split("\n")
    new_text = ""
    new_text += wrapper.fill(text_lines[0]).rstrip("\n") + "\n"
    wrapper.initial_indent = " " * leading_space
    for text_line in text_lines[1:]:
        new_text += wrapper.fill(text_line).rstrip("\n") + "\n"
    return new_text


def parse(val):
    t = html_parser.handle(val).strip() if isinstance(val, str) else val
    return t

# def parse(val):
#     ss = re.sub(r'<br>', '\n', val)
#     ss = re.sub(r'<.*?>', '', ss)
#     return ss

def parse_dict(obj):
    if isinstance(obj, dict):
        for key in obj:
            if isinstance(obj[key], str):
                obj[key] = parse(obj[key])
            elif isinstance(obj[key], dict):
                parse_dict(obj[key])
    return obj

def gs(dic, *args, **kwargs):
    if 'default' not in kwargs:
        kwargs['default'] = None

    d = dic
    for a in args:
        if a in d:
            d = d[a]
        else:
            return kwargs['default']
    if isinstance(d, str) or isinstance(d, unicode):
        d = strip_html(d)
    elif isinstance(d, list):
        for i, element in enumerate(d):
            if isinstance(element, str) or isinstance(element, unicode):
                d[i] = strip_html(element)
    return d

def strip_html(text):
    pat = re.compile("<html>(.+?)</html>", re.S)
    text = pat.sub("\g<1>", text)
    pat = re.compile("<body>(.+?)</body>", re.S)
    text = pat.sub("\g<1>", text)
    pat = re.compile("<div>(.+?)</div>", re.S)
    text = pat.sub("\g<1>\n", text)
    pat = re.compile("<b>(.+?)</b>", re.S)
    text = pat.sub("\g<1>", text)
    pat = re.compile("<i>(.+?)</i>", re.S)
    text = pat.sub("\g<1>", text)
    text = re.sub("<br>", "\n", text)
    text = re.sub("&nbsp;", " ", text)
    text = re.sub("\n{2,}", "\n\n", text)
    return text

def snapgene_file_to_dict(filepath=None, fileobject=None):
    """Return a dictionnary containing the data from a ``*.dna`` file.

    Parameters
    ----------
    filepath
      Path to a .dna file created with SnapGene
    fileobject
      On object-like pointing to the data of a .dna file created with SnapGene
    """

    if filepath is not None:
        fileobject = open(filepath, 'rb')

    if fileobject.read(1) != b'\t':
        raise ValueError('Wrong format for a SnapGene file !')

    def unpack(size, mode):
        return struct.unpack('>' + mode, fileobject.read(size))[0]

    # READ THE DOCUMENT PROPERTIES
    length = unpack(4, 'I')
    title = fileobject.read(8).decode('ascii')
    if length != 14 or title != 'SnapGene':
        raise ValueError('Wrong format for a SnapGene file !')

    data = dict(
        isDNA=unpack(2, 'H'),
        exportVersion=unpack(2, 'H'),
        importVersion=unpack(2, 'H'),
        features=[]
    )

    while True:
        # READ THE WHOLE FILE, BLOCK BY BLOCK, UNTIL THE END
        next_byte = fileobject.read(1)

        # next_byte table
        # 0: dna sequence
        # 1: compressed DNA
        # 2: unknown
        # 3: unknown
        # 5: primers
        # 6: notes
        # 7: history tree
        # 8: additional sequence properties segment
        # 9: file Description
        # 10: features
        # 11: history node
        # 13: unknown
        # 16: alignable sequence
        # 17: alignable sequence
        # 18: sequence trace
        # 19: Uracil Positions
        # 20: custom DNA colors

        if next_byte == b'':
            # END OF FILE
            break

        block_size = unpack(4, 'I')

        if ord(next_byte) == 0:
            # READ THE SEQUENCE AND ITS PROPERTIES
            props = unpack(1, 'b')
            data["dna"] = dict(
                topology="circular" if props & 0x01 else "linear",
                strandedness="double" if props & 0x02 > 0 else "single",
                damMethylated=props & 0x04 > 0,
                dcmMethylated=props & 0x08 > 0,
                ecoKIMethylated=props & 0x10 > 0,
                length=block_size - 1
            )
            data["seq"] = fileobject.read(block_size - 1).decode('ascii')

        elif ord(next_byte) == 6:
            # READ THE NOTES
            block_content = fileobject.read(block_size).decode('utf-8')
            note_data = parse_dict(xmltodict.parse(block_content))
            data['notes'] = note_data['Notes']

        elif ord(next_byte) == 10:
            # READ THE FEATURES
            strand_dict = {"0": ".", "1": "+", "2": "-", "3": "="}
            format_dict = {'@text': parse, '@int': int}
            features_data = xmltodict.parse(fileobject.read(block_size))

            for f in features_data["Features"]["Feature"]:
                segments = f["Segment"]
                if not isinstance(segments, list):
                    segments = [segments]
                segments_ranges = [
                    sorted([int(e) for e in segment['@range'].split('-')])
                    for segment in segments
                ]
                qualifiers = f.get('Q', [])
                if not isinstance(qualifiers, list):
                    qualifiers = [qualifiers]
                parsed_qualifiers = {}
                for q in qualifiers:
                    if isinstance(q['V'], list):
                        if len(q['V'][0].items()) == 1:
                            parsed_qualifiers[q['@name']] = l = []
                            for e in q['V']:
                                fmt, value = e.popitem()
                                fmt = format_dict.get(fmt, parse)
                                l.append(fmt(value))
                        else:
                            parsed_qualifiers[q['@name']] = d = {}
                            for e in q['V']:
                                (fmt1, value1), (fmt2, value2) = e.items()
                                fmt = format_dict.get(fmt1, parse)
                                d[value2] = fmt(value1)

                    else:
                        fmt, value = q['V'].popitem()
                        fmt = format_dict.get(fmt, parse)
                        parsed_qualifiers[q['@name']] = fmt(value)

                if 'label' not in parsed_qualifiers:
                    parsed_qualifiers['label'] = f['@name']
                if 'note' not in parsed_qualifiers:
                    parsed_qualifiers['note'] = []
                if not isinstance(parsed_qualifiers['note'], list):
                    parsed_qualifiers['note'] = [parsed_qualifiers['note']]
                color = segments[0]['@color']
                parsed_qualifiers['note'].append("color: " + color)

                data["features"].append(dict(
                    start=min([start for (start, end) in segments_ranges]),
                    end=max([end for (start, end) in segments_ranges]),
                    strand=strand_dict[f.get('@directionality', "0")],
                    type=f['@type'],
                    name=f['@name'],
                    color=segments[0]['@color'],
                    textColor='black',
                    segments=segments,
                    row=0,
                    isOrf=False,
                    qualifiers=parsed_qualifiers
                ))

        else:
            # WE IGNORE THE WHOLE BLOCK
            # fileobject.read(block_size)
            # 2:
            # 3:

            block = fileobject.read(block_size)
            print(ord(next_byte), len(block))
            # data['block_{}'.format(ord(next_byte))] = block

    fileobject.close()

    return data

def snapgene_file_to_seqrecord(filepath=None, fileobject=None):
    """Return a BioPython SeqRecord from the data of a ``*.dna`` file.

    Parameters
    ----------
    filepath
        Path to a .dna file created with SnapGene
    fileobject
        On object-like pointing to the data of a .dna file created with SnapGene
    """
    data = snapgene_file_to_dict(filepath=filepath, fileobject=fileobject)
    strand_dict = {'+': 1, '-': -1, '.': 0}

    return SeqRecord(
        seq=Seq(data['seq'], alphabet=DNAAlphabet()),
        features=[
            SeqFeature(
                location=FeatureLocation(
                    start=feature['start'],
                    end=feature['end'],
                    strand=strand_dict[feature['strand']]
                ),
                strand=strand_dict[feature['strand']],
                type=feature['type'],
                qualifiers= feature['qualifiers']
            )
            for feature in data['features']
        ],
        annotations=dict(data['notes'])
    )

def snapgene_file_to_gbk(read_file_object, write_file_object):
    data = snapgene_file_to_dict(fileobject=read_file_object)
    with open('dst.json', 'w') as jf:
        jf.write(json.dumps(data, indent=4))
    w = write_file_object
    if data["isDNA"]:
        sequence_unit = "bp"
        sequence_type = "DNA"
    else:
        sequence_unit = "aa"
        sequence_type = ""
    w.write('LOCUS       Exported                {} {} {}s-{}     {} SYN {}\n'\
            .format(gs(data, 'dna', 'length'), sequence_unit,
                    data['dna']['strandedness'][0], sequence_type,
                    data['dna']['topology'], time.strftime("%d-%b-%Y", time.localtime()).upper()))
    definition = gs(data, 'notes', 'Description', default='.').replace('\n', '\n            ')
    w.write(wrap_text('DEFINITION  {}\n'.format(definition), 12))
    w.write('ACCESSION   .\n')
    w.write('VERSION     .\n')
    w.write(wrap_text('KEYWORDS    {}\n'.format(gs(data, 'notes', 'CustomMapLabel', default='.')), 12))
    w.write('SOURCE      .\n')
    w.write('  ORGANISM  .\n')

    references = gs(data, 'notes', 'References')

    reference_count = 0
    if references:
        for key in references:
            reference_count += 1
            ref = references[key]
            w.write('REFERENCE   {}  (bases 1 to {} )\n'.format(reference_count, gs(data, 'dna', 'length')))
            for key2 in ref:
                gb_key = key2.replace('@', '').upper()
                w.write('  {:10}{}\n'.format(gb_key, strip_html(ref[key2])))

    # generate special reference
    reference_count += 1
    w.write('REFERENCE   {}  (bases 1 to {} )\n'.format(reference_count, gs(data, 'dna', 'length')))
    w.write('  AUTHORS   SnapGeneReader\n')
    w.write('  TITLE     Direct Submission\n')
    w.write('  JOURNAL   Exported from SnapGene File Reader\n')
    w.write(wrap_text('COMMENT     {}\n'.format(gs(data, 'notes', 'Comments', default='.'). \
                       replace('\\', '')), 12))
    w.write('FEATURES             Location/Qualifiers\n')

    features = gs(data, 'features')
    features = sorted(features, key=lambda x: gs(x, "start", default=9999999))
    for feature in features:
        strand = gs(feature, 'strand', default='')
        segments = gs(feature, 'segments', default=[])
        segments = [x for x in segments if x['@type'] == 'standard']
        if len(segments) > 1:
            line = 'join('
            seg_positions = []
            for segment in segments:
                segment_range = gs(segment, '@range').replace('-','..')
                seg_match = re.search(r"^(\d+)\.\.(\d+)$", segment_range)
                if seg_match:
                    seg_positions.append((int(seg_match.group(1)), int(seg_match.group(2))))
                if gs(segment, '@type') == 'standard':
                    line += segment_range
                    line += ','
            line = line[:-1] + ')'
            countinous_flag = True
            for i in range(len(seg_positions) - 1):
                if seg_positions[i][1] + 1 != seg_positions[i+1][0]:
                    countinous_flag = False
            if countinous_flag:
                line = '{}..{}'.format(seg_positions[0][0], seg_positions[-1][1])

        else:
            line = '{}..{}'.format(
                gs(feature, 'start', default=' '),
                gs(feature, 'end', default=' ')
            )

        if strand == '-':
            w.write(wrap_text('     {} complement({})\n'.format(
                gs(feature, 'type', default=' ').ljust(15),
                line,
                ), 21))
        else:
            w.write(wrap_text('     {} {}\n'.format(
                gs(feature, 'type', default=' ').ljust(15),
                line,
                ), 21))
        strand = gs(feature, 'strand', default='')
        # if strand == '-':
        #     w.write('                     /direction=LEFT\n')
        # name
        w.write(wrap_text('                     /label="{}"\n'.format(
            gs(feature, 'name', default='.')
        ), 21))
        written_name = gs(feature, 'name', default='.')
        # qualifiers
        for q_key in gs(feature, 'qualifiers', default={}):
            # do not write label, because it has been written at first.
            # if q_key == 'label':
            #     pass
            if q_key == 'note':
                for note in gs(feature, 'qualifiers', q_key, default=[]):
                    # do note write color, because it will be written later
                    if note[:6] != 'color:':
                        w.write(wrap_text('                     /note="{}"\n'.format(note), 21))
            elif q_key == 'direction':
                # do not write direction
                pass
            elif q_key == 'label':
                # write the label if there are multiple 'label' qualifiers
                labels = gs(feature, 'qualifiers', q_key, default=[])
                if isinstance(labels, list):
                    for label in labels:
                        if label != written_name:
                            #do not write name again
                            w.write(wrap_text('                     /{}="{}"\n'.format(
                                    q_key, label), 21))
                else:
                    pass
            else:
                quals = gs(feature, 'qualifiers', q_key, default='')
                if isinstance(quals, list):
                    for qual in quals:
                        if q_key == "translation":
                            qual = re.sub(r",", "", qual)
                        w.write(wrap_text('                     /{}="{}"\n'.format(
                                q_key, qual), 21))
                else:
                    if q_key == "translation":
                        quals = re.sub(r",", "", quals)
                    w.write(wrap_text('                     /{}="{}"\n'.format(
                            q_key, quals), 21))
        # if len(segments) > 1:
        #     w.write('                     /note="This feature has {} segments:'.format(len(segments)))
        #     for seg_i in range(len(segments)):
        #         segment_name = gs(segments[seg_i], '@name', default='')
        #         if segment_name:
        #             segment_name = ' / {}'.format(segment_name)
        #         w.write('\n                        {}:  {} / {}{}'.format(
        #                 seg_i,
        #                 segments[seg_i]['@range'].replace('-',' .. '),
        #                 segments[seg_i]['@color'],
        #                 segment_name,
        #             )
        #         )
        #     w.write('"\n')
        # else:
        #     pass
        #     write colors and direction
        #     w.write('                     /note="color: {}'.format(gs(feature, 'color', default='#ffffff')))
        #     if strand == '-':
        #         w.write('; direction: LEFT"\n')
        #         # w.write('"\n')
        #     elif strand == '+':
        #         w.write('; direction: RIGHT"\n')
        #     else:
        #         w.write('"\n')

    # sequence
    w.write('ORIGIN\n')
    seq = gs(data, 'seq')
    # devide rows
    for i in range(0, len(seq), 60):
        w.write(str(i+1).rjust(9))
        for j in range(i, min(i+60, len(seq)), 10):
            w.write(' {}'.format(seq[j:j+10]))
        w.write('\n')
    w.write('//\n')


# with open('tests/test_samples/custom.dna', 'rb') as f:
#     with open('dst.gbk', 'w', encoding='utf-8') as w:
#         snapgene_file_to_gbk(f, w)
