# Function to translate IDs into results
def ids_to_vals(ids):
    if isinstance(ids, ray.ObjectID):
        ids = ray.get(ids)
    if isinstance(ids, ray.ObjectID):
        return ids_to_vals(ids)
    if isinstance(ids, list):
        results = []
        for id in ids:
            results.append(ids_to_vals(id))
        return results
    return ids