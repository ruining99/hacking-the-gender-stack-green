from rest_framework.serializers import Serializer
from rest_framework.serializers import CharField


class ValidRequestSerializer(Serializer):
    smiles = CharField(required=True)
