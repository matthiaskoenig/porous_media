"""
This file FROM InterFEBIO https://github.com/andresutrera/interFEBio
from enum import IntEnum, auto, Enum
Modified for VTK by Jacob Sturdy 2023.
"""
from enum import Enum, IntEnum, auto


class FEDataType(IntEnum):
    FLOAT = 0  # // scalar             : single fp
    VEC3F = auto()  # // 3D vector          : 3 fps
    MAT3FS = auto()  # // symm 2o tensor     : 6 fps
    MAT3FD = auto()  # // diagonal 2o tensor : 3 fps
    TENS4FS = auto()  # // symm 4o tensor     : 21 fps
    MAT3F = auto()  # // 2o tensor          : 9 fps
    # PLT_ARRAY		= auto() #// variable array (see dictionary for size)
    # PLT_ARRAY_VEC3F	= auto() #// array of vec3f (see dictionary for size)


class FEDataDim(IntEnum):
    FLOAT = 1
    VEC3F = 3
    MAT3FS = 6
    MAT3FD = 3
    TENS4FS = 21
    MAT3F = 9
    # PLT_ARRAY		=
    # PLT_ARRAY_VEC3F	=


class Storage_Fmt(IntEnum):
    FMT_NODE = 0
    FMT_ITEM = auto()
    FMT_MULT = auto()
    FMT_REGION = auto()
    FMT_MATPOINTS = auto()


class Elem_Type(IntEnum):
    ELEM_HEX = 0
    ELEM_PENTA = auto()
    ELEM_TET4 = auto()
    ELEM_QUAD = auto()
    ELEM_TRI = auto()
    ELEM_LINE2 = auto()
    ELEM_HEX20 = auto()
    ELEM_TET10 = auto()
    ELEM_TET15 = auto()
    ELEM_HEX27 = auto()
    ELEM_TRI6 = auto()
    ELEM_QUAD8 = auto()
    ELEM_QUAD9 = auto()
    ELEM_PENTA15 = auto()
    ELEM_TET20 = auto()
    ELEM_TRI10 = auto()
    ELEM_PYRA5 = auto()
    ELEM_TET5 = auto()
    ELEM_PYRA13 = auto()


class nodesPerElementClass(IntEnum):
    ELEM_HEX = 8
    ELEM_PENTA = 6
    ELEM_TET4 = 4
    ELEM_QUAD = 4
    ELEM_TRI = 3
    ELEM_LINE2 = 2
    ELEM_HEX20 = 20
    ELEM_TET10 = 10
    ELEM_TET15 = 15
    ELEM_HEX27 = 27
    ELEM_TRI6 = 6
    ELEM_QUAD8 = 8
    ELEM_QUAD9 = 9
    ELEM_PENTA15 = 15
    ELEM_TET20 = 20
    ELEM_TRI10 = 10
    ELEM_PYRA5 = 5
    ELEM_TET5 = 5
    ELEM_PYRA13 = 13


class tags(Enum):
    PLT_VERSION = "0x0031"
    PLT_ROOT = "0x01000000"
    PLT_HEADER = "0x01010000"
    PLT_HDR_VERSION = "0x01010001"
    # //	PLT_HDR_NODES = '0x01010002' ,
    # //	PLT_HDR_MAX_FACET_NODES = '0x01010003' ,	#// removed (redefined in seach SURFACE section)
    PLT_HDR_COMPRESSION = "0x01010004"
    PLT_HDR_AUTHOR = "0x01010005"  # 	// new in 2.0
    PLT_HDR_SOFTWARE = "0x01010006"  # // new in 2.0
    PLT_DICTIONARY = "0x01020000"
    PLT_DIC_ITEM = "0x01020001"
    PLT_DIC_ITEM_TYPE = "0x01020002"
    PLT_DIC_ITEM_FMT = "0x01020003"
    PLT_DIC_ITEM_NAME = "0x01020004"
    PLT_DIC_ITEM_ARRAYSIZE = "0x01020005"  # // added in version 0x05
    PLT_DIC_ITEM_ARRAYNAME = "0x01020006"  # // added in version 0x05
    PLT_DIC_GLOBAL = "0x01021000"
    # //	PLT_DIC_MATERIAL	 = '0x01022000' ,#	// this was removed
    PLT_DIC_NODAL = "0x01023000"
    PLT_DIC_DOMAIN = "0x01024000"
    PLT_DIC_SURFACE = "0x01025000"
    # //PLT_MATERIALS	 = '0x01030000' ,	#	// This was removed
    # //	PLT_MATERIAL = '0x01030001' ,
    ##//	PLT_MAT_ID	 = '0x01030002' ,
    # //	PLT_MAT_NAME = '0x01030003' ,
    PLT_MESH = "0x01040000"  # 	// this was PLT_GEOMETRY
    PLT_NODE_SECTION = "0x01041000"
    PLT_NODE_HEADER = "0x01041100"  # 	// new in 2.0
    PLT_NODE_SIZE = "0x01041101"  # 	// new in 2.0
    PLT_NODE_DIM = "0x01041102"  # 	// new in 2.0
    PLT_NODE_NAME = "0x01041103"  # 	// new in 2.0
    PLT_NODE_COORDS = "0x01041200"  # 	// new in 2.0
    PLT_DOMAIN_SECTION = "0x01042000"
    PLT_DOMAIN = "0x01042100"
    PLT_DOMAIN_HDR = "0x01042101"
    PLT_DOM_ELEM_TYPE = "0x01042102"
    PLT_DOM_PART_ID = "0x01042103"  # // this was PLT_DOM_MAT_ID
    PLT_DOM_ELEMS = "0x01032104"
    PLT_DOM_NAME = "0x01032105"
    PLT_DOM_ELEM_LIST = "0x01042200"
    PLT_ELEMENT = "0x01042201"
    PLT_SURFACE_SECTION = "0x01043000"
    PLT_SURFACE = "0x01043100"
    PLT_SURFACE_HDR = "0x01043101"
    PLT_SURFACE_ID = "0x01043102"
    PLT_SURFACE_FACES = "0x01043103"
    PLT_SURFACE_NAME = "0x01043104"
    PLT_SURFACE_MAX_FACET_NODES = (
        "0x01043105"  # // new in 2.0 (max number of nodes per facet)
    )
    PLT_FACE_LIST = "0x01043200"
    PLT_FACE = "0x01043201"
    PLT_NODESET_SECTION = "0x01044000"
    PLT_NODESET = "0x01044100"
    PLT_NODESET_HDR = "0x01044101"
    PLT_NODESET_ID = "0x01044102"
    PLT_NODESET_NAME = "0x01044103"
    PLT_NODESET_SIZE = "0x01044104"
    PLT_NODESET_LIST = "0x01044200"
    PLT_PARTS_SECTION = "0x01045000"  # // new in 2.0
    PLT_PART = "0x01045100"
    PLT_PART_ID = "0x01045101"
    PLT_PART_NAME = "0x01045102"
    # 	// plot objects were added in 3.0
    PLT_OBJECTS_SECTION = "0x01050000"
    PLT_OBJECT_ID = "0x01050001"
    PLT_OBJECT_NAME = "0x01050002"
    PLT_OBJECT_TAG = "0x01050003"
    PLT_OBJECT_POS = "0x01050004"
    PLT_OBJECT_ROT = "0x01050005"
    PLT_OBJECT_DATA = "0x01050006"
    PLT_POINT_OBJECT = "0x01051000"
    PLT_POINT_COORD = "0x01051001"
    PLT_LINE_OBJECT = "0x01052000"
    PLT_LINE_COORDS = "0x01052001"
    PLT_STATE = "0x02000000"
    PLT_STATE_HEADER = "0x02010000"
    PLT_STATE_HDR_ID = "0x02010001"
    PLT_STATE_HDR_TIME = "0x02010002"
    PLT_STATE_STATUS = "0x02010003"  # // new in 3.1
    PLT_STATE_DATA = "0x02020000"
    PLT_STATE_VARIABLE = "0x02020001"
    PLT_STATE_VAR_ID = "0x02020002"
    PLT_STATE_VAR_DATA = "0x02020003"
    PLT_GLOBAL_DATA = "0x02020100"
    # //PLT_MATERIAL_DATA = '0x02020200' ,// this was removed
    PLT_NODE_DATA = "0x02020300"
    PLT_ELEMENT_DATA = "0x02020400"
    PLT_FACE_DATA = "0x02020500"
    PLT_MESH_STATE = "0x02030000"
    PLT_ELEMENT_STATE = "0x02030001"
    PLT_OBJECTS_STATE = "0x02040000"
