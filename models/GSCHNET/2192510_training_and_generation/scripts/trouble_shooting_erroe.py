import ase.db
import json

def fix_metadata(db_path):
    db = ase.db.connect(db_path)
    try:
        # Extract existing metadata
        metadata = json.loads(db.metadata.get('metadata', '{}'))
        
        # Check if the property units are stored within metadata
        if '_property_unit_dict' in metadata:
            # Extract the correct dictionary and set it as a top-level key
            db.metadata['_property_unit_dict'] = json.dumps(metadata['_property_unit_dict'])
            print("Successfully moved _property_unit_dict to top-level metadata.")
        else:
            print("No '_property_unit_dict' found in the existing metadata.")
    except Exception as e:
        print("Error during metadata extraction:", e)

db_path = "/home/chem/msummc/GeoLDM_drug/GeoLDM/data/geom/geom_drugs_30.db"
fix_metadata(db_path)




# def check_metadata(db_path):
#     db = ase.db.connect(db_path)
#     try:
#         metadata = db.metadata
#         print("Metadata:", metadata)
#         if '_property_unit_dict' in metadata:
#             print("Property Units Found:", json.loads(metadata['_property_unit_dict']))
#         else:
#             print("No '_property_unit_dict' found.")
#     except Exception as e:
#         print("Error reading metadata:", e)

# db_path = "/home/chem/msummc/GeoLDM_drug/GeoLDM/data/geom/geom_drugs_30.db"
# check_metadata(db_path)
