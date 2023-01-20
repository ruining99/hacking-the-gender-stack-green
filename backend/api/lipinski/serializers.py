from rest_framework.serializers import Serializer
from rest_framework.serializers import CharField


class LipinskiRequestSerializer(Serializer):
    smiles = CharField(required=True)
