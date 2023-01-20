from rest_framework.response import Response
from rest_framework.views import APIView
from api.number.serializers import NumberRequestSerializer
from drf_spectacular.utils import OpenApiExample
from drf_spectacular.utils import extend_schema
from api.number.examples import EXAMPLE_REQUEST
from science.rdkit_endpoints import valid


class Number(APIView):
    serializer_class = NumberRequestSerializer

    @extend_schema(
        examples=[
            OpenApiExample(name="Example", value=EXAMPLE_REQUEST, request_only=True)
        ]
    )
    def post(self, request):
        serializer = NumberRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        smiles = serializer.validated_data["smiles"]
        num = valid(smiles)
        print("science returned,", num)
        print("hi")
        return Response(num)
