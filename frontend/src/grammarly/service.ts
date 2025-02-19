import ApiError from '../shared/xhr/ApiError';
import fetch from '../shared/xhr/fetch';

async function checkValidity(smiles: string) {
  const response = await fetch('api/valid/', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

async function checkProps(smiles: string) {
  const response = await fetch('api/number/', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

async function checkLipinski(smiles: string) {
  const response = await fetch('api/lipinski/', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

async function getImage(smiles: string) {
  const response = await fetch('api/image/', {
    method: 'GET',
    body: JSON.stringify({ smiles }),
  });

  if (!response.ok) {
    throw new ApiError(response.clone());
  }

  return response.json();
}

declare global {
  interface Window {
    checkValidity: any;
    getImage: any;
    checkProps: any;
    checkLipinski: any;
    submit: any;
    props: any;
    lipinski: any;
  }
}
window.checkValidity = checkValidity;
window.checkLipinski = checkLipinski;
window.checkProps = checkProps;
window.getImage = getImage;
export { checkLipinski, checkValidity, checkProps, getImage };

// async function registerNewRGroup(smiles: string) {
//   const response = await fetch(rgroupApiUrls.create, {
//     method: 'POST',
//     body: JSON.stringify({ smiles }),
//   });

//   if (!response.ok) {
//     throw new ApiError(response.clone());
//   }

//   return response.json() as Promise<RGroup>;
// }

// async function enumerateProperties(enumerationOpts: {
//   core_smiles: string;
//   rgroup_smiles: Record<string, string[]>;
// }) {
//   const response = await fetch(enumerationApiUrls.create, {
//     method: 'POST',
//     body: JSON.stringify(enumerationOpts),
//   });

//   if (!response.ok) {
//     throw new ApiError(response.clone());
//   }

//   return response.json() as Promise<EnumerateResponse>;
// }

// export type { RGroup, EnumerateResponse };
// export { fetchAllRGroups, registerNewRGroup, enumerateProperties };
