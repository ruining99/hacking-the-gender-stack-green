from rest_framework.response import Response
from rest_framework.views import APIView
from api.valid.serializers import ValidRequestSerializer
from drf_spectacular.utils import OpenApiExample
from drf_spectacular.utils import extend_schema
from api.valid.examples import EXAMPLE_REQUEST
from science.rdkit_endpoints import valid


class Valid(APIView):
    serializer_class = ValidRequestSerializer

    @extend_schema(
        examples=[
            OpenApiExample(name="Example", value=EXAMPLE_REQUEST, request_only=True)
        ]
    )
    def post(self, request):
        serializer = ValidRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        smiles = serializer.validated_data["smiles"]
        isvalid = valid(smiles)
        print("science returned,", isvalid)
        print("hi")
        return Response(isvalid)
