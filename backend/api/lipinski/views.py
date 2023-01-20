from rest_framework.response import Response
from rest_framework.views import APIView
from api.lipinski.serializers import LipinskiRequestSerializer
from drf_spectacular.utils import OpenApiExample
from drf_spectacular.utils import extend_schema
from api.lipinski.examples import EXAMPLE_REQUEST
from science.rdkit_endpoints import lipinski_violations
from science.rdkit_endpoints import lipinski_test_result


class Lipinski(APIView):
    serializer_class = LipinskiRequestSerializer

    @extend_schema(
        examples=[
            OpenApiExample(name="Example", value=EXAMPLE_REQUEST, request_only=True)
        ]
    )
    def post(self, request):
        serializer = LipinskiRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        smiles = serializer.validated_data["smiles"]

        # valid function returns a boolean
        lipinski = lipinski_test_result(smiles)
        print("science returned,", lipinski)
        return Response(lipinski)
