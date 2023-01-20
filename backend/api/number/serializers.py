from rest_framework.serializers import Serializer
from rest_framework.serializers import CharField


class NumberRequestSerializer(Serializer):
    smiles = CharField(required=True)
